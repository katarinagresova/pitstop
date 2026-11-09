"""Read repair (`pitstop repair`): clips the untemplated 5' nt. It writes a new BAM in which the extra 5' nucleotides that
do not match the transcript are soft clips. The other commands skip a 5' soft clip (bam.py), so they read the repaired
BAM with no change.

The clip looks at the first UNTEMPLATED_MAX nt of a read (inspect.UNTEMPLATED_MAX). From the first aligned nt on, a nt
that does not match the transcript becomes a soft clip; it stops at the first match, at an N (in the read or the
transcript), and keeps at least one nt aligned. A 5' soft clip of the aligner counts towards the UNTEMPLATED_MAX nt, so
a second run clips nothing. `inspect` uses the same rule for its values after repair.

repair looks at the alignments the other commands read (bam.py: forward, a CIGAR of matches and soft clips at the ends)
on GTF transcripts, at each alignment of a read. It needs the read sequence: of the alignment, or of the alignment of
the same read just before it (a secondary alignment can have none). It changes the position and the CIGAR, and NM and
MD where the alignment has them; not the sequence, the qualities or other tags (AS stays). Other records are copied.

Outputs:
  <prefix>.bam         in the order of the input. A BAM sorted by coordinate stays sorted (a clip moves an alignment
                       by up to UNTEMPLATED_MAX nt) and gets an index, <prefix>.bam.bai. The header gets a @PG line.
  <prefix>.repair.tsv  per footprint length (before the repair: the read length minus the 5' soft clip, as `inspect`)
                       and clipped nt (0..UNTEMPLATED_MAX): the alignments that repair looked at, and the reads (their
                       primary alignments)
"""
import heapq
import logging
import os
import re
from collections import Counter

import pandas as pd
import pysam

from . import __version__, annotation, bam
from .inspect import UNTEMPLATED_MAX
from .tables import write_tsv

log = logging.getLogger(__name__)

ACGT = frozenset("ACGT")
_MD = re.compile(r"\d+|\^[A-Za-z]+|[A-Za-z]")


def shape(cigar):
    """(footprint length, 5' soft clip, aligned nt) of a CIGAR that the other commands read, else None."""
    fp = bam._footprint_length(cigar)
    if not fp:
        return None
    c5 = cigar[0][1] if cigar[0][0] == bam._SOFT else 0
    return fp, c5, sum(n for op, n in cigar if op in bam._MATCH)


def to_clip(seq, c5, aligned, tx, pos):
    """The untemplated 5' nt to clip of a read (seq) whose first aligned nt (seq[c5]) is at tx[pos]."""
    k = 0
    while c5 + k < UNTEMPLATED_MAX and k < aligned - 1:
        r, t = seq[c5 + k], tx[pos + k]
        if r == t or r not in ACGT or t not in ACGT:
            break
        k += 1
    return k


def clip_cigar(cigar, k):
    """The CIGAR with its first k aligned nt soft-clipped (k below the aligned nt)."""
    c5 = cigar[0][1] if cigar[0][0] == bam._SOFT else 0
    rest = list(cigar[1 if c5 else 0:])
    left = k
    while left:
        op, n = rest[0]
        if n > left:
            rest[0], left = (op, n - left), 0
        else:
            rest.pop(0)
            left -= n
    return [(bam._SOFT, c5 + k)] + rest


def clip_md(md, k):
    """(the MD tag without its first k aligned nt, the mismatches among them). The first k nt have no deletion."""
    tok = _MD.findall(md)
    i, left, mism = 0, k, 0
    while left:
        n = int(tok[i])
        if n >= left:
            tok[i] = str(n - left)
            break
        if tok[i + 1].startswith("^"):
            raise ValueError(f"MD {md}: a deletion in the first {k} nt")
        left -= n + 1
        mism += 1
        i += 2
    return "".join(tok[i:]), mism


def clip(a, k):
    """Soft-clip the first k aligned nt of alignment a (in place)."""
    a.cigartuples = clip_cigar(a.cigartuples, k)
    a.reference_start += k
    tags = a.get_tags(with_value_type=True)
    names = [t[0] for t in tags]
    if "MD" not in names and "NM" not in names:
        return
    mism = k
    if "MD" in names:
        j = names.index("MD")
        md, mism = clip_md(tags[j][1], k)
        tags[j] = ("MD", md, "Z")
    if "NM" in names:
        j = names.index("NM")
        tags[j] = ("NM", max(tags[j][1] - mism, 0), tags[j][2])
    a.set_tags(tags)


def pg_line(header):
    """The header dict with a @PG line for this run (after the last one)."""
    pg = header.setdefault("PG", [])
    ids = {p.get("ID") for p in pg}
    pid, n = "pitstop", 0
    while pid in ids:
        n += 1
        pid = f"pitstop.{n}"
    line = {"ID": pid, "PN": "pitstop", "VN": __version__,
            "DS": f"repair: up to {UNTEMPLATED_MAX} untemplated 5' nt (no match to the transcript) become a soft clip"}
    if pg:
        line["PP"] = pg[-1]["ID"]
    pg.append(line)
    return header


def load_transcripts(fa, fasta, refs, gtf_tx, ref_len):
    """Per BAM reference: its transcript sequence, None if not in the GTF. Loaded in genome order: in the order of the
    alignments of an unsorted BAM, the FASTA was read at random places, and that took most of the time."""
    rank = {c: j for j, c in enumerate(fa.references)}
    for r, t in zip(refs, gtf_tx):
        if t is not None and t["chrom"] not in rank:
            raise ValueError(f"{t['chrom']} (transcript {r}) is not in {fasta}")
    seqs = [None] * len(refs)
    for i in sorted((i for i, t in enumerate(gtf_tx) if t is not None),
                    key=lambda i: (rank[gtf_tx[i]["chrom"]], min(s for s, _ in gtf_tx[i]["exons"]))):
        seqs[i] = annotation._tx_seq(fa, gtf_tx[i])
        if len(seqs[i]) != ref_len[i]:
            raise ValueError(f"{refs[i]} is {ref_len[i]} nt in the BAM but {len(seqs[i])} nt in the GTF")
    return seqs


def run(bam_path, gtf, fasta, out_prefix):
    out_bam = f"{out_prefix}.bam"
    if os.path.abspath(out_bam) == os.path.abspath(bam_path):
        raise ValueError(f"{out_bam} is the input BAM: give another --out-prefix")
    txs = annotation.read_gtf(gtf)
    counts = Counter()        # (footprint length, clipped nt, primary) -> alignments
    skipped = Counter()
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as inp:
        refs = [r.split("|", 1)[0] for r in inp.references]
        gtf_tx = [txs.get(r) for r in refs]
        if not any(t is not None for t in gtf_tx):
            raise ValueError(f"no BAM reference is a GTF transcript (BAM: {', '.join(refs[:3])}; "
                             f"GTF: {', '.join(list(txs)[:3])}): do the transcript ids differ?")
        seqs = load_transcripts(annotation._open_fasta(fasta), fasta, refs, gtf_tx, inp.lengths)
        shapes = {}           # CIGAR string -> shape(), () if not read by the other commands
        header = pg_line(inp.header.to_dict())
        by_coord = header.get("HD", {}).get("SO") == "coordinate"
        buffer, tid, n = [], None, 0   # by_coord: (new position, input order, alignment) of the reference tid
        last_name = last_seq = None
        with pysam.AlignmentFile(out_bam, "wb", header=header) as out:
            for a in inp.fetch(until_eof=True):
                pos = a.reference_start
                k = None
                if a.flag & bam._SKIP:
                    skipped["unmapped, reverse or supplementary"] += 1
                else:
                    c = a.cigarstring
                    sh = shapes.get(c)
                    if sh is None:
                        sh = shapes[c] = shape(a.cigartuples) or ()
                    i = a.reference_id
                    if not sh:
                        skipped["CIGAR"] += 1
                    elif gtf_tx[i] is None:
                        skipped["reference not in the GTF"] += 1
                    else:
                        fp, c5, aligned = sh
                        seq = a.query_sequence
                        name = a.query_name
                        if seq is None and name == last_name:
                            seq = last_seq
                        elif seq is not None:
                            last_name, last_seq = name, seq
                        if seq is None:
                            skipped["no read sequence"] += 1
                        else:
                            k = to_clip(seq, c5, aligned, seqs[i], pos)
                            counts[fp, k, not a.flag & 0x100] += 1
                            if k:
                                clip(a, k)
                if not by_coord:
                    out.write(a)
                    continue
                if a.reference_id != tid:
                    for x in sorted(buffer):
                        out.write(x[2])
                    buffer, tid = [], a.reference_id
                while buffer and buffer[0][0] <= pos:      # later alignments start at pos or after
                    out.write(heapq.heappop(buffer)[2])
                heapq.heappush(buffer, (a.reference_start, n, a))
                n += 1
            for x in sorted(buffer):
                out.write(x[2])
    if by_coord:
        pysam.index(out_bam)
    table = pd.DataFrame([(fp, k, counts[fp, k, True] + counts[fp, k, False], counts[fp, k, True])
                          for fp, k in sorted({(fp, k) for fp, k, _ in counts})],
                         columns=["length", "clipped", "alignments", "reads"])
    write_tsv(table, f"{out_prefix}.repair.tsv")
    seen = int(table["alignments"].sum())
    changed = table.loc[table["clipped"] > 0, ["alignments", "reads"]].sum()
    log.info("%s: repair clipped %d of %d alignments (%.2f%%), %d reads; not looked at: %s", bam_path,
             changed["alignments"], seen, 100 * changed["alignments"] / max(seen, 1), changed["reads"],
             ", ".join(f"{v} {k}" for k, v in sorted(skipped.items())) or "none")
    return table

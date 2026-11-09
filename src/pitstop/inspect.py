"""Read QC of one library (`pitstop inspect`): which read lengths keep one reading frame, where their reads fall, what
is at the read ends, how complex the library is; flags for common library problems, and the read lengths to use.

The reads are those of the other commands (bam.read_bam): forward alignments whose CIGAR is matches plus soft clips
at the ends; a 5' soft clip is not part of the footprint. Alignments to BAM references that are not GTF transcripts
are counted and dropped (as in `orfs`). A read is split evenly over its alignments to transcripts with a CDS or, if it
has none, over all its alignments (as quant's qc.tsv), so each read counts once. A length is a footprint length: the
read length minus the 5' soft clip.

Words. The CDS interior: the read is at least INTERIOR nt inside the CDS at both ends. The phase of a read: (5' end -
CDS start) mod 3; the top frame of a length: the share of its top phase in the CDS interior. The start peak: the
commonest 5' end at the annotated starts (nt from the first nt of the start codon), in PEAK; beside: the larger of the
two nt beside a peak, as a share of the peak. The offset of a length: minus its start peak, or DEFAULT_OFFSET with
fewer than MIN_ANCHOR reads at the peak; its P-sites give the regions (5' UTR, CDS, 3' UTR: the stop codon and after;
transcript: no CDS). Repair (`pitstop repair`) clips the untemplated 5' nt: of the first UNTEMPLATED_MAX nt of the
read, from the first aligned nt on, each that does not match the transcript becomes a soft clip; it stops at the
first match (an aligner's 5' soft clip counts towards the UNTEMPLATED_MAX nt, and one nt stays aligned). The 3' end
spread: the standard deviation of the length of the reads whose 5' end is at the start peak of all lengths above
SHORT_MAX and below DISOME (the 5' end is fixed, so the length gives the 3' end). A nuclease that trims every
footprint to the ribosome gives about one length there; under-digestion leaves the 3' ends at many places.

A class per length, the first rule that holds:
  few reads    fewer than MIN_FRAME reads in the CDS interior, or less than MIN_SHARE of the reads
  short state  SHORT nt, the start peak at SHORT_PEAK with MIN_ANCHOR reads or more, SHORT_CDS of the P-sites or more in
               the CDS: ribosomes with an open A site (not damaged footprints; repair does not touch them)
  clean        top frame CLEAN or more
  no frame     else
The class uses the frame in the CDS interior, as quant, orfs and score do. Mixed (a column apart from the class): a clean
or no-frame length with MIN_ANCHOR reads or more at the start peak, and beside it MIXED of them or more (as `starts`):
two groups of 5' ends at the starts, and one offset puts the smaller group 1 nt off. It matters for the start lengths
of `starts` only.
The suggested read lengths: of the runs of consecutive clean lengths, the one with the most reads.

Tables (<prefix>.<name>.tsv). The report (report/inspect.py) reads only these, so every number of the report is in a table.
  lengths     per length: the reads and their share; the CDS-interior reads and their phases; the start peak (its reads,
              share and beside); the 3' ends at the stops the same way (nt from the first nt of the stop codon); the
              offset; the P-sites by region (shares); the class and mixed; after repair, the reads, top frame and
              class
  anchors     per length: the reads by 5' end around the annotated starts (START_WINDOW), and by 3' end around the
              stops (STOP_WINDOW); share: of the length's reads in that window
  ends        per length and end, and for all lengths (length "all"): the reads, their soft clips (1, 2, 3 or more nt),
              the mismatch rate at positions 0..K-1 from that end (the soft clip included: position 0 is the first nt
              of the read), the read nt at a mismatched position 0, the transcript nt at the cut (the first or the
              last nt of the footprint) and its bias: its share over the share of the nt one codon inside the
              footprint (same codon position, not at a cut)
  complexity  per subsample (crc32 of the read name, mod 10000, below fraction x 10000): the reads and the distinct
              footprints (BAM reference, 5' end, length; of a multimapper, its smallest alignment). The distinct
              share falls with the depth, so the flag uses summary's distinct_share_at_depth: the expected distinct
              footprints per read in a subsample of COMPLEXITY_DEPTH reads (exact, from the reads per footprint; NA
              with fewer reads)
  summary     one row: the inputs, the read counts, the library values (with the 3' end spread: end3_peak, the start
              peak; end3_reads, its reads; end3_spread), the flags (FLAGS) and the suggestions
The flags have one of five states: pass, warn, fail, unknown (too few reads or no read sequences), info (no cut).
"""
import logging
import os
import zlib
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pysam

from . import __version__, annotation, bam
from .annotation import CDS, INTERIOR, NO_CDS, REGIONS, UTR3, UTR5
# as in starts: a class needs MIN_FRAME CDS-interior reads; the start peak needs MIN_ANCHOR reads; mixed: beside the
# start peak MIXED of its reads or more
from .starts import MIN_ANCHOR, MIN_FRAME, MIXED, peak_beside
from .tables import write_tsv

log = logging.getLogger(__name__)

K = bam.K            # positions from each read end in ends.tsv
MIN_SHARE = 0.003    # a class also needs this share of the reads: a few positions can give a tiny length any frame
CLEAN = 0.7          # clean: top frame this or more
DEFAULT_OFFSET = 12  # the offset of a length with fewer than MIN_ANCHOR reads at the start peak
SHORT = (20, 22)     # the short state: these lengths,
SHORT_PEAK = -12     # this start peak,
SHORT_CDS = 0.5      # and this share of the P-sites in the CDS or more
PEAK = (-20, -4)     # the start peak is in this window, and leaves 2 codons (P and A site) in the read
START_WINDOW = (-30, 0)   # anchors.tsv: the 5' ends, nt from the first nt of the start codon
STOP_WINDOW = (-5, 30)    # anchors.tsv: the 3' ends, nt from the first nt of the stop codon
UNTEMPLATED_MAX = 2  # repair clips at most this many untemplated 5' nt
LONG = 32            # long reads: this many nt or more
SHORT_MAX = 22       # short reads: this many nt or fewer
DISOME = (50, 70)    # disomes
BACKGROUND = (4, 8)  # the mismatch background: positions 4-7 from the end
SUBSAMPLES = 10      # complexity: fractions 0.1, 0.2, ... 1
CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
BASES = "ACGT"

# The cuts of the flags, set on 20 public libraries of known problems and 81 clean ones (and their yeast spike-ins).
WEAK_FRAME = (0.5, 0.3)      # clean share: pass >= 0.5, warn >= 0.3, else fail
LONG_READS = (0.10, 0.25)    # share at LONG nt or more: pass < 0.10, warn < 0.25, else fail
UNDER_DIGESTION = (1.4, 1.7)  # 3' end spread (nt): pass < 1.4, warn < 1.7, else fail
SPREAD_READS = 1_000         # ... unknown with fewer reads at the start peak (a GTF without 5' UTRs has none)
UNTEMPLATED = (0.05, 5)      # 5' mismatch >= 0.05 and >= 5x its background ...
CLIPPED = 0.5                # ... and clipped by the aligner (info, not fail) if the 5' clip share >= 0.5 x it
SPIKE = (3, 0.1)             # a mismatch spike: > 3x the background and > 0.1 (5' positions 2-3, 3' positions 0-1)
SHORT_READS = 10_000         # the short state: unknown with fewer reads at SHORT nt
COMPLEXITY_DEPTH = 1_000_000  # low complexity: the distinct share in a subsample of this many reads (unknown with fewer)
COMPLEXITY = 0.4              # ... warn below it. No fail: without UMIs, PCR copies and ribosomes on one nt look alike
FLAGS = ("weak_frame", "long_reads", "under_digestion", "untemplated", "end_spike", "mixed_lengths", "short_state",
         "low_complexity")


@dataclass
class Context:
    """Per alignment: the weight, the CDS (cs, ce: -1 without), the BAM reference; the transcript sequences."""
    w: np.ndarray
    coding: np.ndarray
    cs: np.ndarray
    ce: np.ndarray
    ref: np.ndarray
    flat: np.ndarray      # codes (CODE) of the transcripts with alignments, end to end
    first: np.ndarray     # per BAM reference: its first position in flat
    size: np.ndarray      # per BAM reference: its length (0 if not loaded)

    def nt(self, pos, sel=slice(None)):
        """The transcript nt code (4: none) at transcript position pos of the alignments sel."""
        ref = self.ref[sel]
        n = self.size[ref]
        ok = (pos >= 0) & (pos < n)
        return np.where(ok, self.flat[self.first[ref] + np.clip(pos, 0, np.maximum(n - 1, 0))], 4).astype(np.uint8)


def read_weights(aln):
    """Per alignment: 1/n over the read's n alignments to transcripts with a CDS, or over all n if it has none; 0 for
    its other alignments."""
    has_cds = np.zeros(aln.n_reads, dtype=bool)
    has_cds[aln.read[aln.tx >= 0]] = True
    use = (aln.tx >= 0) | ~has_cds[aln.read]
    n = np.bincount(aln.read[use], minlength=aln.n_reads)
    return np.where(use, 1.0 / np.maximum(n[aln.read], 1), 0.0)


def transcripts(aln, txs, fasta, bam_len):
    """flat, first, size (Context) for the BAM references with alignments."""
    used = np.unique(aln.ref)
    fa = annotation._open_fasta(fasta)
    first = np.zeros(len(aln.refs), dtype=np.int64)
    size = np.zeros(len(aln.refs), dtype=np.int64)
    parts, at = [], 0
    for i in used:
        s = annotation._tx_seq(fa, txs[aln.refs[i]])
        if len(s) != bam_len[i]:
            raise ValueError(f"{aln.refs[i]} is {bam_len[i]} nt in the BAM but {len(s)} nt in the GTF")
        first[i], size[i] = at, len(s)
        parts.append(s)
        at += len(s)
    flat = CODE[np.frombuffer("".join(parts).encode(), dtype=np.uint8)]
    return flat, first, size


def profile(x, li, mask, w, window, n_len):
    """Weighted reads per (length index, x in window), mask: alignments that count."""
    k = mask & (x >= window[0]) & (x <= window[1])
    n = window[1] - window[0] + 1
    return np.bincount(li[k] * n + (x[k] - window[0]), weights=w[k], minlength=n_len * n).reshape(n_len, n)


def peaks(prof, window, allowed):
    """Per row: the position of the peak (the first if tied), its reads, and the larger nt beside it as a share of the
    peak. allowed: (rows, positions) bool."""
    i, at, before, after = peak_beside(prof, allowed)
    beside = np.maximum(before, after)
    has = at > 0
    return (np.where(has, i + window[0], 0), at, np.where(has, beside / np.where(has, at, 1), np.nan), has)


def length_table(pos5, length, ctx, lengths):
    """lengths.tsv's columns for one view of the alignments (pos5, length: int64 per alignment)."""
    lo, hi = lengths
    n_len = hi - lo + 1
    li = length - lo
    w, coding, cs, ce = ctx.w, ctx.coding, ctx.cs, ctx.ce
    p3 = pos5 + length - 1
    L = np.arange(lo, hi + 1)
    reads = np.bincount(li, weights=w, minlength=n_len)
    inner = coding & (pos5 >= cs + INTERIOR) & (p3 < ce - INTERIOR)
    phases = np.bincount(li[inner] * 3 + (pos5 - cs)[inner] % 3, weights=w[inner], minlength=n_len * 3).reshape(-1, 3)
    interior = phases.sum(axis=1)
    start = profile(pos5 - cs, li, coding, w, START_WINDOW, n_len)
    stop = profile(p3 - ce, li, coding, w, STOP_WINDOW, n_len)
    x = np.arange(START_WINDOW[0], START_WINDOW[1] + 1)
    allowed = (x[None, :] >= np.maximum(PEAK[0], -(L[:, None] - 6))) & (x[None, :] <= PEAK[1])
    s_peak, s_at, s_beside, s_has = peaks(start, START_WINDOW, allowed)
    t_peak, t_at, t_beside, t_has = peaks(stop, STOP_WINDOW, np.ones(stop.shape, dtype=bool))
    offset = np.where(s_at >= MIN_ANCHOR, -s_peak, DEFAULT_OFFSET)
    psite = pos5 + offset[li]
    region = np.where(coding, np.select([psite < cs, psite <= ce - 3], [UTR5, CDS], UTR3), NO_CDS)
    n_reg = len(REGIONS)
    regions = np.bincount(li * n_reg + region, weights=w, minlength=n_len * n_reg).reshape(-1, n_reg)
    share = lambda a, b: np.where(b > 0, a / np.where(b > 0, b, 1), np.nan)
    reg = share(regions, reads[:, None])
    top = share(phases.max(axis=1), interior)
    s_start = start.sum(axis=1)
    t_stop = stop.sum(axis=1)
    cls = np.select(
        [(interior < MIN_FRAME) | (reads < MIN_SHARE * reads.sum()),
         (L >= SHORT[0]) & (L <= SHORT[1]) & (s_peak == SHORT_PEAK) & (s_at >= MIN_ANCHOR)
         & (np.nan_to_num(reg[:, CDS]) >= SHORT_CDS),
         np.nan_to_num(top) >= CLEAN],
        ["few reads", "short state", "clean"], "no frame")
    mixed = np.isin(cls, ["clean", "no frame"]) & (s_at >= MIN_ANCHOR) & (np.nan_to_num(s_beside) >= MIXED)
    out = pd.DataFrame({"length": L, "reads": reads, "share": reads / max(reads.sum(), 1e-300), "interior": interior})
    for f in range(3):
        out[f"phase{f}"] = share(phases[:, f], interior)
    out["top_phase"] = pd.array(np.where(interior > 0, phases.argmax(axis=1), 0), dtype="Int64")
    out.loc[interior == 0, "top_phase"] = pd.NA
    out["top_frame"] = top
    for name, total, pk, at, bes, has in (("start", s_start, s_peak, s_at, s_beside, s_has),
                                         ("stop", t_stop, t_peak, t_at, t_beside, t_has)):
        out[f"{name}_reads"] = total
        out[f"{name}_peak"] = pd.array(pk, dtype="Int64")
        out.loc[~has, f"{name}_peak"] = pd.NA
        out[f"{name}_peak_share"] = share(at, total)
        out[f"{name}_beside"] = bes
    out["offset"] = offset
    for i, r in enumerate(REGIONS):
        out[r] = reg[:, i]
    out["class"] = cls
    out["mixed"] = mixed
    return out, start, stop


def untemplated_clip(aln, ends, ctx):
    """Per alignment: the untemplated 5' nt that repair clips (0..UNTEMPLATED_MAX; repair.to_clip for one alignment)."""
    c5 = ends.clip5.astype(np.int64)
    aligned = aln.length.astype(np.int64) - ends.clip3
    k = np.zeros(aln.read.size, dtype=np.int64)
    go = np.ones(aln.read.size, dtype=bool)
    for j in range(UNTEMPLATED_MAX):
        q = c5 + j
        rb = CODE[ends.head[aln.read, np.minimum(q, K - 1)]]
        tb = ctx.nt(aln.pos5.astype(np.int64) + j)
        go &= (q < UNTEMPLATED_MAX) & (j < aligned - 1) & (rb < 4) & (tb < 4) & (rb != tb)
        k += go
    return k


def ends_table(aln, ends, ctx, lengths):
    """ends.tsv: per length and end, and for all lengths (first). Counts are weighted reads."""
    lo, hi = lengths
    n_len = hi - lo + 1
    sel = np.flatnonzero(ctx.w > 0)
    w = ctx.w[sel]
    L = aln.length[sel].astype(np.int64)
    li = L - lo
    read = aln.read[sel]
    p5 = aln.pos5[sel].astype(np.int64)
    rows = []
    for end in (5, 3):
        clip = (ends.clip5 if end == 5 else ends.clip3)[sel].astype(np.int64)
        c5 = ends.clip5[sel].astype(np.int64)
        qlen = c5 + L                      # the read length
        raw0 = p5 - c5                     # the transcript position of the read's first nt
        clips = np.bincount(li * 4 + np.minimum(clip, 3), weights=w, minlength=n_len * 4).reshape(-1, 4)
        n = np.zeros((n_len, K))
        x = np.zeros((n_len, K))
        mm_base = np.zeros((n_len, 4))
        for i in range(K):
            if end == 5:
                rb, tpos = CODE[ends.head[read, i]], raw0 + i
            else:
                rb, tpos = CODE[ends.tail[read, K - 1 - i]], raw0 + qlen - 1 - i
            tb = ctx.nt(tpos, sel)
            ok = (rb < 4) & (tb < 4) & (i < qlen)
            n[:, i] = np.bincount(li[ok], weights=w[ok], minlength=n_len)
            bad = ok & (rb != tb)
            x[:, i] = np.bincount(li[bad], weights=w[bad], minlength=n_len)
            if i == 0:
                mm_base = np.bincount(li[bad] * 4 + rb[bad], weights=w[bad], minlength=n_len * 4).reshape(-1, 4)
        cut_pos, inner_pos = (p5, p5 + 3) if end == 5 else (p5 + L - 1, p5 + L - 4)
        cut, inner = ctx.nt(cut_pos, sel), ctx.nt(inner_pos, sel)
        ok = (cut < 4) & (inner < 4)
        cuts = np.bincount(li[ok] * 4 + cut[ok], weights=w[ok], minlength=n_len * 4).reshape(-1, 4)
        inners = np.bincount(li[ok] * 4 + inner[ok], weights=w[ok], minlength=n_len * 4).reshape(-1, 4)
        reads = np.bincount(li, weights=w, minlength=n_len)
        for name, sl in (("all", slice(None)),) + tuple((str(lo + j), slice(j, j + 1)) for j in range(n_len)):
            r = reads[sl].sum()
            if r == 0:
                continue
            c = clips[sl].sum(axis=0)
            nn, xx = n[sl].sum(axis=0), x[sl].sum(axis=0)
            mb = mm_base[sl].sum(axis=0)
            cu, inn = cuts[sl].sum(axis=0), inners[sl].sum(axis=0)
            cu_s = cu / cu.sum() if cu.sum() else np.full(4, np.nan)
            in_s = inn / inn.sum() if inn.sum() else np.full(4, np.nan)
            row = {"length": name, "end": end, "reads": r, "clip1": c[1] / r, "clip2": c[2] / r, "clip3plus": c[3] / r}
            row.update({f"mm{i}": xx[i] / nn[i] if nn[i] else np.nan for i in range(K)})
            row.update({f"mm0_{b}": mb[j] / mb.sum() if mb.sum() else np.nan for j, b in enumerate(BASES)})
            row.update({f"cut_{b}": cu_s[j] for j, b in enumerate(BASES)})
            row.update({f"bias_{b}": cu_s[j] / in_s[j] if in_s[j] > 0 else np.nan for j, b in enumerate(BASES)})
            rows.append(row)
    out = pd.DataFrame(rows)
    order = np.argsort([(-1 if r == "all" else int(r)) * 2 + (e == 3) for r, e in zip(out["length"], out["end"])],
                       kind="stable")
    return out.iloc[order].reset_index(drop=True)


def anchors_table(start, stop, lengths):
    """anchors.tsv: long format, every length and position."""
    lo, hi = lengths
    parts = []
    for anchor, end, arr, window in (("start", 5, start, START_WINDOW), ("stop", 3, stop, STOP_WINDOW)):
        n = window[1] - window[0] + 1
        total = arr.sum(axis=1, keepdims=True)
        share = np.where(total > 0, arr / np.where(total > 0, total, 1), np.nan)
        parts.append(pd.DataFrame({"length": np.repeat(np.arange(lo, hi + 1), n), "anchor": anchor, "end": end,
                                   "x": np.tile(np.arange(window[0], window[1] + 1), hi - lo + 1),
                                   "reads": arr.ravel(), "share": share.ravel()}))
    return pd.concat(parts, ignore_index=True)


def complexity_table(aln):
    """complexity.tsv: per subsample, the reads and the distinct footprints; and the reads per footprint (all reads).
    The footprint of a read is its smallest (BAM reference, 5' end, length), so the copies of a read (same sequence,
    same alignments) have one footprint."""
    pos = aln.pos5.astype(np.int64)
    length = aln.length.astype(np.int64)
    key = (aln.ref.astype(np.int64) * (int(pos.max()) + 1) + pos) * (int(length.max()) + 1) + length
    o = np.lexsort((key, aln.read))
    r = aln.read[o]
    first = np.r_[True, r[1:] != r[:-1]]
    ids, read_key = r[first], key[o][first]
    names = aln.read_names
    u = np.fromiter((zlib.crc32(names[i].encode()) % 10000 for i in ids), dtype=np.int64, count=ids.size)
    o = np.lexsort((u, read_key))
    k = read_key[o]
    head = np.r_[True, k[1:] != k[:-1]]
    first_u = u[o][head]              # per footprint: the smallest u of its reads
    rows = []
    for j in range(1, SUBSAMPLES + 1):
        cut = j * 10000 // SUBSAMPLES
        rows.append((j / SUBSAMPLES, int((u < cut).sum()), int((first_u < cut).sum())))
    reads_per_footprint = np.diff(np.r_[np.flatnonzero(head), k.size])
    return pd.DataFrame(rows, columns=["fraction", "reads", "distinct"]), reads_per_footprint


def expected_distinct(counts, n):
    """The expected distinct footprints in a subsample of n of the N reads (without replacement, exact), from the reads
    per footprint: the sum over footprints of 1 - P(none of its c reads drawn), P = prod_{j<c} (N - j - n) / (N - j).
    NaN if n > N."""
    N = int(counts.sum())
    if n > N:
        return np.nan
    c, m = np.unique(counts, return_counts=True)
    j = np.arange(min(int(c.max()), N - n))         # a footprint with more than N - n reads is always drawn
    log_none = np.r_[0.0, np.cumsum(np.log1p(-n / (N - j)))]
    drawn = np.where(c <= N - n, -np.expm1(log_none[np.minimum(c, j.size)]), 1.0)
    return float((m * drawn).sum())


def runs(lengths):
    """Runs of consecutive lengths."""
    out = []
    for x in lengths:
        if out and out[-1][-1] == x - 1:
            out[-1].append(x)
        else:
            out.append([x])
    return out


def suggest(table, cls="class", reads="reads"):
    """The suggested read lengths: of the runs of consecutive clean lengths, the one with the most reads (the first if
    tied), as LO-HI; "" if none."""
    by = table.set_index("length")[reads]
    rs = runs(table.loc[table[cls] == "clean", "length"].tolist())
    if not rs:
        return ""
    best = max(rs, key=lambda run: by[run].sum())
    return f"{best[0]}-{best[-1]}"


def flags(s):
    """The state of each flag (FLAGS) from the summary values s."""
    f = {}
    if s["classed_lengths"] == 0:
        f["weak_frame"] = "unknown"
    else:
        c = s["clean_share"]
        f["weak_frame"] = "pass" if c >= WEAK_FRAME[0] else "warn" if c >= WEAK_FRAME[1] else "fail"
    v = s["long_share"]
    f["long_reads"] = "pass" if v < LONG_READS[0] else "warn" if v < LONG_READS[1] else "fail"
    v = s["end3_spread"]
    f["under_digestion"] = ("unknown" if pd.isna(v) else "pass" if v < UNDER_DIGESTION[0] else
                            "warn" if v < UNDER_DIGESTION[1] else "fail")
    m, bg = s["mismatch5"], s["mismatch5_background"]
    if pd.isna(m) or pd.isna(bg):
        f["untemplated"] = "unknown"
    elif m >= UNTEMPLATED[0] and m >= UNTEMPLATED[1] * bg:
        f["untemplated"] = "info" if s["clip5_share"] >= CLIPPED * m else "fail"
    else:
        f["untemplated"] = "pass"
    spikes = [(s["mismatch5_inner"], s["mismatch5_background"]), (s["mismatch3"], s["mismatch3_background"])]
    if any(pd.isna(a) or pd.isna(b) for a, b in spikes):
        f["end_spike"] = "unknown"
    else:
        f["end_spike"] = "warn" if any(a > SPIKE[0] * b and a > SPIKE[1] for a, b in spikes) else "pass"
    f["mixed_lengths"] = "info" if s["mixed_lengths"] else "pass"
    if s["short_reads"] < SHORT_READS:
        f["short_state"] = "unknown"
    else:
        f["short_state"] = "warn" if s["short_cds_share"] < SHORT_CDS else "info"
    d = s["distinct_share_at_depth"]
    f["low_complexity"] = "unknown" if pd.isna(d) else "warn" if d < COMPLEXITY else "pass"
    return f


def end3_spread(anchors):
    """(start peak, its reads, 3' end spread) of the lengths above SHORT_MAX and below DISOME, from anchors.tsv; the
    spread is NaN with fewer than SPREAD_READS reads at the peak."""
    a = anchors[(anchors["anchor"] == "start") & (anchors["length"] > SHORT_MAX) & (anchors["length"] < DISOME[0])
                & (anchors["x"] >= PEAK[0]) & (anchors["x"] <= PEAK[1])]
    by_x = a.groupby("x")["reads"].sum()
    if by_x.empty or by_x.max() <= 0:
        return np.nan, 0.0, np.nan
    peak = int(by_x.idxmax())
    at = a[a["x"] == peak]
    w, L = at["reads"].to_numpy(), at["length"].to_numpy()
    n = w.sum()
    if n < SPREAD_READS:
        return peak, n, np.nan
    m = (w * L).sum() / n
    return peak, n, np.sqrt((w * (L - m) ** 2).sum() / n)


def tv_distance(a, b):
    """Total-variation distance of two profiles; NaN if one is empty."""
    if a.sum() <= 0 or b.sum() <= 0:
        return np.nan
    return 0.5 * np.abs(a / a.sum() - b / b.sum()).sum()


def summary(stats, aln, lengths_df, anchors_df, ends_df, cplx, reads_per_footprint, names):
    """summary.tsv's one row (a dict)."""
    s = dict(names)
    s.update({k: stats[k] for k in ("reads_forward", "alignments_dropped_cigar", "reads_cigar_ok", "refs_not_in_gtf",
                                    "alignments_dropped_ref_not_in_gtf")})
    nh = np.bincount(aln.read, minlength=aln.n_reads)
    nh = nh[nh > 0]
    s["reads"] = int(nh.size)
    on_cds = np.zeros(aln.n_reads, dtype=bool)
    on_cds[aln.read[aln.tx >= 0]] = True
    s["reads_on_cds_transcript"] = int(on_cds.sum())
    s.update(reads_nh1=int((nh == 1).sum()), reads_nh2=int((nh == 2).sum()), reads_nh3=int((nh == 3).sum()),
             reads_nh4plus=int((nh >= 4).sum()))
    full, before = cplx.iloc[-1], cplx.iloc[-2]
    s["distinct_share"] = full["distinct"] / full["reads"] if full["reads"] else np.nan
    s["new_per_read"] = ((full["distinct"] - before["distinct"]) / (full["reads"] - before["reads"])
                         if full["reads"] > before["reads"] else np.nan)
    s["distinct_share_at_depth"] = expected_distinct(reads_per_footprint, COMPLEXITY_DEPTH) / COMPLEXITY_DEPTH
    t = lengths_df
    L, share = t["length"].to_numpy(), t["share"].to_numpy()
    s["common_length"] = int(L[np.argmax(t["reads"].to_numpy())])
    s["max_length"] = int(L[t["reads"].to_numpy() > 0].max())
    s["long_share"] = share[L >= LONG].sum()
    s["short_share"] = share[L <= SHORT_MAX].sum()
    s["disome_share"] = share[(L >= DISOME[0]) & (L <= DISOME[1])].sum()
    s["end3_peak"], s["end3_reads"], s["end3_spread"] = end3_spread(anchors_df)
    reads = t["reads"].to_numpy()
    cds = np.nan_to_num(t["CDS"].to_numpy()) * reads
    no_cds = np.nan_to_num(t["transcript"].to_numpy()) * reads
    s["cds_share"] = cds.sum() / reads.sum()
    sh = (L >= SHORT[0]) & (L <= SHORT[1])
    s["short_reads"] = reads[sh].sum()
    s["short_cds_share"] = cds[sh].sum() / reads[sh].sum() if reads[sh].sum() else np.nan
    s["length_distance"] = tv_distance(cds, no_cds)
    s["classed_lengths"] = int((t["class"] != "few reads").sum())
    s["clean_share"] = share[t["class"] == "clean"].sum()
    s["clean_lengths"] = ",".join(str(x) for x in L[t["class"] == "clean"])
    s["read_lengths"] = suggest(t)
    rep = t["reads_repaired"].to_numpy()
    s["clean_share_repaired"] = rep[t["class_repaired"] == "clean"].sum() / max(rep.sum(), 1e-300)
    s["clean_lengths_repaired"] = ",".join(str(x) for x in L[t["class_repaired"] == "clean"])
    s["read_lengths_repaired"] = suggest(t, "class_repaired", "reads_repaired")
    framed = t[~t["class"].isin(["few reads", "short state"])]
    if len(framed):
        b = framed.loc[framed["top_frame"].idxmax()]
        s["best_length"], s["best_top_frame"] = int(b["length"]), b["top_frame"]
    else:
        s["best_length"], s["best_top_frame"] = pd.NA, np.nan
    s["mixed_lengths"] = ",".join(str(x) for x in L[t["mixed"].to_numpy(dtype=bool)])
    a5 = ends_df[(ends_df["length"] == "all") & (ends_df["end"] == 5)].iloc[0]
    a3 = ends_df[(ends_df["length"] == "all") & (ends_df["end"] == 3)].iloc[0]
    bg = [f"mm{i}" for i in range(*BACKGROUND)]
    s["mismatch5"] = a5["mm0"]
    s["mismatch5_inner"] = np.fmax(a5["mm2"], a5["mm3"])
    s["mismatch5_background"] = a5[bg].astype(float).mean()
    s["mismatch5_T"] = a5["mm0_T"]
    s["mismatch3"] = np.fmax(a3["mm0"], a3["mm1"])
    s["mismatch3_background"] = a3[bg].astype(float).mean()
    s["clip5_share"] = a5["clip1"] + a5["clip2"] + a5["clip3plus"]
    # an added nt that matches the transcript is not seen: about 1/4 of the transcript nt are T, so the hidden part is
    # about 1/3 of the seen part
    s["untemplated_hidden"] = max(0.0, s["mismatch5"] - s["mismatch5_background"]) / 3 if pd.notna(s["mismatch5"]) \
        else np.nan
    bias = [(abs(np.log2(r[f"bias_{b}"])), r[f"bias_{b}"], f"{e}' {b}")
            for e, r in ((5, a5), (3, a3)) for b in BASES if r[f"bias_{b}"] > 0]
    top = max(bias, key=lambda x: x[0]) if bias else (np.nan, np.nan, "")
    s["bias_max"], s["bias_max_at"] = top[1], top[2]
    f = flags(s)
    s.update({f"flag_{k}": f[k] for k in FLAGS})
    s["repair"] = "yes" if f["untemplated"] == "fail" else "no"
    return s


def write(df, path):
    write_tsv(df.round(6), path, na_rep="NA")


def run(bam_path, gtf, fasta, out_prefix, report=True):
    stats = {}
    txs = annotation.read_gtf(gtf)
    anno = annotation.load_annotation(gtf, fasta, txs=txs)
    aln, ends = bam.read_bam(bam_path, anno, stats, keep=set(txs), with_ends=True)
    if aln.read.size == 0:
        raise ValueError(f"{bam_path}: no alignment to a GTF transcript with an acceptable CIGAR")
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as f:
        bam_len = np.array(f.lengths, dtype=np.int64)
    flat, first, size = transcripts(aln, txs, fasta, bam_len)
    coding = aln.tx >= 0
    t = np.maximum(aln.tx, 0)
    ctx = Context(read_weights(aln), coding, np.where(coding, anno.cds_start[t], -1),
                  np.where(coding, anno.cds_end[t], -1), aln.ref, flat, first, size)
    pos5, length = aln.pos5.astype(np.int64), aln.length.astype(np.int64)
    k = untemplated_clip(aln, ends, ctx)
    lengths = (int(min(length.min(), (length - k).min())), int(length.max()))
    now, start, stop = length_table(pos5, length, ctx, lengths)
    rep, _, _ = length_table(pos5 + k, length - k, ctx, lengths)
    now["reads_repaired"] = rep["reads"]
    now["top_frame_repaired"] = rep["top_frame"]
    now["class_repaired"] = rep["class"]
    keep = (now["reads"] > 0) | (now["reads_repaired"] > 0)
    now = now[keep].reset_index(drop=True)
    anchors = anchors_table(start, stop, lengths)
    anchors = anchors[anchors["length"].isin(now["length"])].reset_index(drop=True)
    ends_df = ends_table(aln, ends, ctx, lengths)
    cplx, reads_per_footprint = complexity_table(aln)
    names = dict(library=os.path.basename(out_prefix), pitstop=__version__, bam=os.path.basename(bam_path),
                 gtf=os.path.basename(gtf), fasta=os.path.basename(fasta))
    s = summary(stats, aln, now, anchors, ends_df, cplx, reads_per_footprint, names)
    write(now, f"{out_prefix}.lengths.tsv")
    write(anchors, f"{out_prefix}.anchors.tsv")
    write(ends_df, f"{out_prefix}.ends.tsv")
    write(cplx, f"{out_prefix}.complexity.tsv")
    write(pd.DataFrame([s]), f"{out_prefix}.summary.tsv")
    log.info("clean share %.3f (after repair: %.3f); read lengths %s; flags: %s", s["clean_share"],
             s["clean_share_repaired"], s["read_lengths"] or "none",
             ", ".join(f"{k} {s['flag_' + k]}" for k in FLAGS if s["flag_" + k] != "pass"))
    if report:
        from . import report as rpt
        rpt.inspect_report(out_prefix)
    return s

"""Footprint alignments from a transcriptome BAM.

Kept: mapped, forward-strand (sense to the transcript), primary or secondary
alignments whose CIGAR is matches plus optional soft clips at the ends. A 5' soft
clip is taken as untemplated nucleotides added during library prep: the 5' end is
the first aligned base and the length is the read length minus the 5' clip.
Alignments with I/D/N/H/P are dropped and counted. Each read name is one read
(names are UMI-collapsed upstream), whatever its number of alignments.
"""
import logging
from array import array
from dataclasses import dataclass

import numpy as np
import pysam

log = logging.getLogger(__name__)

_MATCH = {0, 7, 8}   # M, =, X
_SOFT = 4
_SKIP = 0x4 | 0x10 | 0x800   # unmapped, reverse, supplementary


@dataclass
class Alignments:
    """The per-alignment arrays are int32 (half the memory of int64). A product of two of them, or of one
    and a large number, can overflow int32: cast to int64 first."""
    read: np.ndarray     # read index (0..n_reads-1), one per alignment
    tx: np.ndarray       # annotation index of the transcript, -1 if it has no valid CDS
    pos5: np.ndarray     # 0-based transcript position of the 5' end
    length: np.ndarray
    n_reads: int
    read_names: list     # read name per read index, in order of first appearance
    ref: np.ndarray = None   # BAM reference index, one per alignment
    refs: list = None        # transcript id per BAM reference

    def subset(self, mask):
        return Alignments(self.read[mask], self.tx[mask], self.pos5[mask], self.length[mask],
                           self.n_reads, self.read_names, self.ref[mask], self.refs)

    def count_reads(self):
        return int(np.unique(self.read).size)


def _footprint_length(cigar):
    """Read length minus the 5' soft clip, or None if the CIGAR is not matches + end soft clips."""
    first, last = 0, len(cigar)
    three = 0
    if cigar[0][0] == _SOFT:
        first = 1
    if last > first and cigar[-1][0] == _SOFT:
        three, last = cigar[-1][1], last - 1
    aligned = 0
    for op, n in cigar[first:last]:
        if op not in _MATCH:
            return None
        aligned += n
    return aligned + three if aligned else None


def read_bam(path, anno, stats, keep=None):
    """All forward alignments with an acceptable CIGAR. stats gets the read counts. With `keep` (transcript
    ids), alignments to other BAM references are counted but not kept: refs_not_in_gtf,
    alignments_dropped_ref_not_in_gtf (those with an acceptable CIGAR). All other counts are as without `keep`."""
    index = anno.index()
    names, last, r = {}, None, -1
    cigar_ok = bytearray()   # per read: 1 if it has an alignment with an acceptable CIGAR
    fp_of = {}               # CIGAR string -> footprint length, 0 if not acceptable
    read, ref, pos5, length = array("i"), array("i"), array("i"), array("i")   # C int: int32
    n_bad_cigar_aln = n_dropped = 0
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        refs = [r.split("|", 1)[0] for r in bam.references]
        ref_tx = np.array([index.get(r, -1) for r in refs], dtype=np.int32)
        ref_len = np.array(bam.lengths, dtype=np.int64)
        known = ref_tx >= 0
        if not known.any():
            raise ValueError(f"no BAM reference is a transcript with a CDS in the GTF (BAM: {', '.join(refs[:3])}; "
                             f"GTF: {', '.join(anno.tx[:3])}): do the transcript ids differ, e.g. in a version suffix?")
        bad = known &(ref_len != anno.tx_len[np.maximum(ref_tx, 0)])
        if bad.any():
            i = np.flatnonzero(bad)[0]
            raise ValueError(f"{bam.references[i]} is {ref_len[i]} nt in the BAM but "
                             f"{anno.tx_len[ref_tx[i]]} nt in the annotation; {bad.sum()} such transcripts")
        stats["cds_transcripts_in_bam_header"] = int(known.sum())
        kept = [keep is None or t in keep for t in refs]
        for a in bam.fetch(until_eof=True):
            if a.flag & _SKIP:
                continue
            q = a.query_name
            if q != last:   # a read's alignments are usually adjacent
                r = names.setdefault(q, len(names))
                if r == len(cigar_ok):
                    cigar_ok.append(0)
                last = q
            c = a.cigarstring
            fp = fp_of.get(c)
            if fp is None:
                fp = fp_of[c] = _footprint_length(a.cigartuples) or 0
            if not fp:
                n_bad_cigar_aln += 1
                continue
            cigar_ok[r] = 1
            i = a.reference_id
            if not kept[i]:
                n_dropped += 1
                continue
            read.append(r)
            ref.append(i)
            pos5.append(a.reference_start)
            length.append(fp)
    stats["reads_forward"] = len(names)
    stats["alignments_dropped_cigar"] = n_bad_cigar_aln
    stats["reads_cigar_ok"] = cigar_ok.count(1)
    if keep is not None:
        stats["refs_not_in_gtf"] = kept.count(False)
        stats["alignments_dropped_ref_not_in_gtf"] = n_dropped
    ref = np.frombuffer(ref, dtype=np.int32)
    aln = Alignments(np.frombuffer(read, dtype=np.int32), ref_tx[ref], np.frombuffer(pos5, dtype=np.int32),
                     np.frombuffer(length, dtype=np.int32), len(names), list(names), ref, refs)
    log.info("%s: %d forward reads, %d with an acceptable alignment (%d alignments dropped for their CIGAR)",
             path, len(names), stats["reads_cigar_ok"], n_bad_cigar_aln)
    return aln

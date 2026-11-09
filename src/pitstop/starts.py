"""Start sites from harringtonine P-sites, against the P-sites of matched elongation libraries.

Harringtonine stops ribosomes at the start codon, and the ribosomes that started before it run off.
A used start has more harringtonine P-sites in its window (the start codon +-1 codon) than the
elongation profile of the same transcript expects. A pause has the same peak in both.

Layout. The transcripts sit end to end in one array, GAP empty nt apart, so that no window reaches
the next transcript. Each position has a region: 5' UTR, CDS, 3' UTR (the stop codon and after) or
transcript (no annotated CDS), and a distance in the region: nt to the CDS start in the 5' UTR, codons
from the CDS start in the CDS, nt from the stop codon in the 3' UTR, nt from the 5' end in a transcript.

Calibration and start lengths. A library can put the P-sites of a read length 1-2 nt off: its nuclease can
leave most 5' ends 1 nt longer than the pool's offsets expect (an initiation library's digestion is not the
elongation library's). Per library and read length: the start anchor is the peak of the P-sites at the
annotated starts (nt -15..32); the frame is the dominant frame of the P-sites well inside the CDSs. The length
is calibrated when its start peak is within MAX_SHIFT nt of nt 0 (with MIN_ANCHOR P-sites or more) and the
frame does not clearly put the P-sites elsewhere (FRAME_GUARD); its P-sites then move by minus the peak. The
frame alone cannot tell a codon (-1 and +2 are both frame 2), the start peak can. A length that is not
calibrated moves by its frame if that is clear, else stays. A length is mixed when a nt beside its start peak has
MIXED of the peak's P-sites or more: its 5' ends split between two offsets, and one offset puts the P-sites of the
other 1 nt off. Which of the two nt is larger changes from library to library. So a length that is mixed in half of
the harringtonine libraries or more is a pair length in every library: its anchor is the pair of the start peak and
the larger nt beside it, and the calibration moves the pair onto nt 0 and 1, the centre of the window; the frame
vetoes it only when it is clearly in neither frame of the pair, and without the anchor a clear frame pair moves onto
frames 0 and 1. The start
lengths come from a QC of each read length in each library: in WINDOW (an 80S footprint is about 28 nt; the window
keeps out short RNA fragments), calibrated in every harringtonine library, and periodic in every
library (MIN_FRAME CDS-interior P-sites or more, MIN_TOP of them in the top frame); of the runs of consecutive such
lengths, the one with the most P-sites. Only their P-sites count, one per alignment. <prefix>.lengths.tsv has the
QC and each length's offset (the orfs run's, plus the shift), <prefix>.profiles.tsv the profiles behind it
(P-sites from the pool's offsets, before the calibration), <prefix>.libraries.tsv the start-length P-sites of each
library after the calibration (by region, by frame inside the CDSs, and the start peak: P-sites at the annotated start
codons per start, over the P-sites per in-frame codon inside the CDSs),
<prefix>.pairs.tsv the read-length distance of each pair (a row of the sheet): half the summed absolute difference
of the two libraries' shares of CDS-interior P-sites per read length (0: the same lengths, 1: none in common; another
digestion at PAIR_DISTANCE or more), <prefix>.start_codons.tsv the start codons of the tests and of the calls per
class. The report (report/starts.py, <prefix>.starts.html) reads only the tables.

Background, per pool. The harringtonine / elongation ratio depends on the region, and in the CDS on
the distance from the start (run-off). h = s_t f_b e: a factor per bin of region and distance (BINS),
fit on null windows, and a scale per transcript, over its background. A null window is the window of
a null codon: no start codon (START_CODONS) within 8 nt, outside the window of an annotated start and
its codons 1-10 (RAISED, never background), with mostly unique reads (processed pseudogenes mirror
their parent's start peak).

The test. A window has h harringtonine and e elongation reads. Its elongation rate has a gamma prior:
mean lam0 (the sum over its nt of the elongation density of their zone in the transcript: 5' UTR, CDS
frame 0, 1 or 2, 3' UTR, transcript), shape a per region (from the null windows). Then h ~ NB(mean R max(e, (a + e) / (a / lam0 + 1)), shape a + e),
with R = f_b s_t: the prior lifts a window with few elongation reads, but does not shrink a pause.
Pools add their means and variances. An extra overdispersion phi per region and bin of
expected reads makes at most 1% of the null windows reach p <= 0.01. q is BH per class: the annotated
starts, and the other candidates per region. A scan calibrates phi twice: after a first call, the null
windows in the codons 1-10 of the called starts (RAISED: run-off, as after an annotated start) are dropped,
and the starts are called again.

Candidates: the starts of an ORF table and the annotated starts, or a scan: each start codon (of the chosen
start codons, default START_CODONS) with harringtonine reads in its window (BH counts every such start codon of
the class: the others have p = 1). Null codons stay clear of all START_CODONS, whichever are chosen.
The codons 1-10 of an annotated start are not called. Then the scan calls starts in the order of p:
a start with q < fdr (FDR) is called if no called start stops it, and it stops the candidates within 8 nt
and in its codons 1-10 (harringtonine raises them). Each called start gives an ORF to its first
in-frame stop, typed against the annotated CDS. The catalogue is the uORFs and uoORFs of the called
starts and every annotated CDS.

The shoulder of a candidate (a flag, not part of the test): log2 harringtonine / elongation P-sites in its codons
1-5, minus the same in the 5 codons before its window. Harringtonine raises the codons 1-5 of a used start; most
peaks at codons that cannot start (null codons) have no shoulder.
"""
import logging
import os
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy.special import betainc, gammainc

from . import __version__, annotation
from .annotation import CDS, INTERIOR, NO_CDS, REGIONS, UTR3, UTR5
from .tables import write_stats, write_tsv

log = logging.getLogger(__name__)

GAP = 64
CLASSES = REGIONS + ("annotated",)   # BH classes: the other candidates by region, the annotated starts
ANNOTATED = len(REGIONS)             # the class of the annotated starts
START_CODONS = ("ATG", "CTG", "GTG", "TTG", "AAG", "ACG", "AGG", "ATA", "ATC", "ATT")
WIN = (-1, 3)        # a start's window: nt from its first nt, the start codon and the nt before it
RAISED = (-3, 33)    # an annotated start's window and its codons 1-10: harringtonine raises them
SHOULDER = (3, 18)   # a start's codons 1-5, and the 5 codons before its window: the shoulder
BEFORE = (-18, -3)
NEAR = (-8, 8)       # a null codon has no start codon this near: their windows would overlap
PROFILE = (-15, 63)  # nt from the annotated start: the kernel to codon 10, then codons 11-20 scale the elongation
KERNEL_END = 33
MAX_SHIFT = 2        # calibration: a start peak this many nt or fewer from nt 0 can be moved onto it
MIN_ANCHOR = 20      # calibration: P-sites at the start peak, summed over the annotated starts
MIXED = 0.5          # a length is mixed when a nt beside its start peak has this share of the peak's P-sites or more
MIN_FRAME = 100      # the frame decides (guard, fallback) only with this many P-sites inside the CDSs
FRAME_GUARD = 1.5    # the frame is clear when its dominant frame has this many times the P-sites of the other
WINDOW = (24, 34)    # start lengths: an 80S footprint is about 28 nt; the window keeps out short RNA fragments
MIN_TOP = 0.4        # start lengths: the top frame has this share of a length's CDS-interior P-sites or more
STOP_PROFILE = (-45, 15)   # QC profiles: nt from the first nt of the stop codon (codons -15..4)
CDS_BINS = 20        # QC profiles: bins along the CDS
LENGTH_COLUMNS = ["library", "pool", "kind", "length", "psites", "share", "5UTR", "CDS", "3UTR", "transcript",
                  "frame0", "frame1", "frame2", "interior_share", "top_frame", "periodic", "start_peak", "at_peak",
                  "beside_peak", "beside_nt", "beside_share", "mixed", "pair", "calibrated", "shift", "offset", "qc",
                  "start_length"]
LIBRARY_COLUMNS = ["library", "pool", "kind", "psites", "start_psites", "start_share", "5UTR", "CDS", "3UTR",
                   "transcript", "interior", "frame0", "frame1", "frame2", "at_start", "start_enrichment"]
PAIR_COLUMNS = ["pool", "harringtonine", "elongation", "distance"]
PAIR_DISTANCE = 0.4  # pairs: a read-length distance this large or larger means another digestion (pairs of one
                     # digestion were at 0.01-0.33, pairs of another digestion at 0.56-0.78)
CODON_COLUMNS = ["class", "codon", "tests", "test_share", "called", "call_share"]
STARTS_COLUMNS = ["Name", "start", "codon", "region", "frame", "ORF_id", "harringtonine", "elongation", "expected",
                  "enrichment", "p", "typical_p", "q", "multimapped", "shoulder"]
SCAN_COLUMNS = ["called", "stopped_by", "end", "type"]   # starts.tsv of a scan: STARTS_COLUMNS, then these
ORF_COLUMNS = ["ORF_id", "Name", "start", "end", "type", "start_codon", "harringtonine", "elongation", "expected",
               "enrichment", "p", "q", "shoulder"]       # start_orfs.tsv, catalogue.tsv
STRIDE = WIN[1] - WIN[0]   # one null codon per window width for the null windows: they do not overlap
FDR = 0.05           # the BH cut of a call
MAX_MULTI = 0.5      # null windows: multimapped share of the window's reads below this
PRIOR_READS = 10     # pseudo-reads: transcript scale at the global scale, bin factor at 1
PRIOR_NT = 100       # pseudo-nt at the region's mean elongation density, per transcript and region
DEPTH = (0, 0.3, 1, 3, 10, 30)               # lower edges of the bins of expected reads in the window
TAIL = 0.01                                  # phi: at most this share of null windows at p <= TAIL
PHI_GRID = np.r_[0, np.geomspace(1e-3, 100, 51)]
MIN_NULL = 200       # null windows for a phi of its own: else the region's, else all windows'
# lower edges of the distance bins per region (the region's distance unit)
BINS = {
    "5UTR": (0, 30, 100, 300),
    "CDS": (0, 20, 50, 100, 200, 400, 700, 1000),
    "3UTR": (0, 100, 300, 1000),
    "transcript": (0,),
}


@dataclass
class Layout:
    names: np.ndarray       # transcript ids, sorted
    first: np.ndarray       # flat position of each transcript's nt 0
    tx_len: np.ndarray
    cds_start: np.ndarray   # -1 without an annotated CDS
    cds_end: np.ndarray     # one past the last sense codon; -1 without an annotated CDS
    seq: np.ndarray         # uint8 ASCII, N in the gaps
    size: int

    def positions(self):
        """Per flat position: transcript index, region index (REGIONS) and distance; -1 in the gaps."""
        tx = np.full(self.size, -1, dtype=np.int32)
        region = np.full(self.size, -1, dtype=np.int8)
        dist = np.full(self.size, -1, dtype=np.int32)
        for i, (f, n, cs, ce) in enumerate(zip(self.first, self.tx_len, self.cds_start, self.cds_end)):
            tx[f:f + n] = i
            x = np.arange(n)
            if cs < 0:
                region[f:f + n], dist[f:f + n] = NO_CDS, x
            else:
                region[f:f + n] = np.select([x < cs, x < ce], [UTR5, CDS], UTR3)
                dist[f:f + n] = np.select([x < cs, x < ce], [cs - x, (x - cs) // 3], x - ce)
        return tx, region, dist

    def cds_starts(self):
        """Flat positions of the annotated starts, and their transcript indices."""
        i = np.flatnonzero(self.cds_start >= 0)
        return self.first[i] + self.cds_start[i], i

    def flat(self, coord):
        """Per transcript: the flat position of coord (a transcript coordinate per transcript, e.g. cds_start); -1
        without an annotated CDS."""
        return np.where(self.cds_start >= 0, self.first + coord, -1)


ZONES = 6   # 5' UTR, CDS frame 0, 1 and 2, 3' UTR, transcript


def zone_index(layout, tx, region):
    """Per flat position: 5' UTR 0, CDS frame 0-2 1-3, 3' UTR 4, transcript 5; -1 in the gaps. The prior
    density of the elongation reads is per nt of a zone: a window holds one or two frame-0 nt of a CDS."""
    z = np.where(region < 0, -1, np.choose(np.maximum(region, 0), [0, 1, 4, 5])).astype(np.int8)
    cds = np.flatnonzero(region == CDS)
    t = tx[cds]
    z[cds] = 1 + (cds - layout.first[t] - layout.cds_start[t]) % 3
    return z


def make_layout(gtf, fasta, names):
    """The layout of the transcripts `names` that are in the GTF."""
    anno = annotation.load_annotation(gtf, fasta)
    seqs = annotation.transcript_seqs(gtf, fasta, names)
    ids = np.array(sorted(seqs), dtype=object)
    idx = anno.index()
    c = np.array([idx.get(t, -1) for t in ids], dtype=np.int64)
    tx_len = np.array([len(seqs[t]) for t in ids], dtype=np.int64)
    first = GAP + np.r_[0, np.cumsum(tx_len + GAP)[:-1]]
    size = int(first[-1] + tx_len[-1] + GAP)
    seq = np.full(size, ord("N"), dtype=np.uint8)
    for t, f in zip(ids, first):
        seq[f:f + len(seqs[t])] = np.frombuffer(seqs[t].encode(), dtype=np.uint8)
    has = c >= 0
    log.info("layout: %d transcripts (%d with a CDS), %d nt", len(ids), has.sum(), size)
    return Layout(ids, first, tx_len, np.where(has, anno.cds_start[c], -1), np.where(has, anno.cds_end[c], -1),
                  seq, size)


def codon_mask(layout, codons=START_CODONS):
    """True at the positions where one of `codons` starts."""
    code = np.full(256, 4, dtype=np.int64)
    for i, b in enumerate(b"ACGT"):
        code[b] = i
    s = code[layout.seq]
    k = 25 * s[:-2] + 5 * s[1:-1] + s[2:]
    table = np.zeros(125, dtype=bool)
    for cod in codons:
        table[25 * "ACGT".index(cod[0]) + 5 * "ACGT".index(cod[1]) + "ACGT".index(cod[2])] = True
    return np.r_[table[k], False, False]


def covered(size, lo, hi):
    """True at the positions in any [lo, hi)."""
    d = np.zeros(size + 1, dtype=np.int64)
    np.add.at(d, lo, 1)
    np.add.at(d, hi, -1)
    return np.cumsum(d)[:-1] > 0


def background(layout, tx):
    """True at the transcript positions outside the RAISED zone of the annotated starts."""
    pos, _ = layout.cds_starts()
    return (tx >= 0) & ~covered(layout.size, pos + RAISED[0], pos + RAISED[1])


def null_codons(bg, starts, near=NEAR):
    """Positions c whose window is all background (bg) and that have no position of `starts` in
    [c + near[0], c + near[1]]. A gap is longer than a window, so the window is in c's transcript."""
    m, s = np.r_[0, np.cumsum(bg)], np.r_[0, np.cumsum(starts)]
    c = np.flatnonzero(bg)
    return c[(m[c + WIN[1]] - m[c + WIN[0]] == WIN[1] - WIN[0]) & (s[c + near[1] + 1] - s[c + near[0]] == 0)]


class PSites(NamedTuple):
    """P-sites, one per alignment."""
    pos: np.ndarray       # flat position
    length: np.ndarray    # read length
    multi: np.ndarray     # the read has another alignment


def read_psites(prefix, layout, lengths):
    """The P-sites in <prefix>.psites.tsv whose length is in [lo, hi] and whose transcript is in the layout."""
    first = pd.Series(layout.first, index=layout.names)
    ps = pd.read_csv(f"{prefix}.psites.tsv", sep="\t", dtype={"read": str, "Name": str, "psite": np.int64,
                                                              "length": np.int64})
    ps = ps[ps["length"].between(*lengths)]
    multi = ps["read"].duplicated(keep=False).to_numpy()
    f = ps["Name"].map(first).to_numpy()
    k = ~np.isnan(f)
    return PSites(f[k].astype(np.int64) + ps["psite"].to_numpy()[k], ps["length"].to_numpy()[k].astype(np.int16),
                  multi[k])


def window_sums(x, pos, win=WIN):
    """Sum of x over [pos + win[0], pos + win[1])."""
    c = np.r_[0, np.cumsum(x)]
    return c[pos + win[1]] - c[pos + win[0]]


def shoulder(pools, pos):
    """Per start at pos: log2 (harringtonine + 1) / (elongation + 1) P-sites in its codons 1-5 (SHOULDER), minus the
    same in the 5 codons before its window (BEFORE). The reads of the pools add, without factors: a local contrast."""
    def log_ratio(win):
        h = sum(window_sums(pl.h, pos, win) for pl in pools)
        e = sum(window_sums(pl.e, pos, win) for pl in pools)
        return np.log2((h + 1) / (e + 1))
    return log_ratio(SHOULDER) - log_ratio(BEFORE)


def bin_index(region, dist):
    """Index of the (region, distance) bin, over BINS in order; -1 in the gaps."""
    out = np.full(region.shape, -1, dtype=np.int32)
    k = 0
    for r, name in enumerate(REGIONS):
        m = region == r
        out[m] = k + np.searchsorted(BINS[name], dist[m], "right") - 1
        k += len(BINS[name])
    return out


def bin_names():
    return [f"{name} {lo}" for name in REGIONS for lo in BINS[name]]


def fit_factors(h, e, tx, b, fit, scale, n_tx, n_bins, prior_reads, iters=500, tol=1e-9):
    """Poisson fit of h = s_t f_b e by alternating updates (h, e, tx, b per position): the factor per bin f
    over the positions `fit`, scaled to an e-weighted mean of 1 there, and the scale per transcript s over
    the positions `scale`. Priors: prior_reads pseudo-reads at the global scale for s; one pseudo-read at
    factor 1 for f, so a bin without reads does not get 0."""
    def sums(x, at):
        return np.bincount(tx[at].astype(np.int64) * n_bins + b[at], weights=x[at],
                           minlength=n_tx * n_bins).reshape(n_tx, n_bins)
    hs, es, hf, ef = sums(h, scale), sums(e, scale), sums(h, fit), sums(e, fit)
    f = np.ones(n_bins)
    for _ in range(iters):
        fe = (es * f).sum(axis=1)
        s = (hs.sum(axis=1) + prior_reads * hs.sum() / fe.sum()) / (fe + prior_reads)
        new = (hf.sum(axis=0) + 1) / ((s[:, None] * ef).sum(axis=0) + 1)
        new /= (new * ef.sum(axis=0)).sum() / ef.sum()
        done = np.abs(new - f).max() < tol
        f = new
        if done:
            break
    return f, s


def nb_sf(h, mean, size):
    """P(X >= h), X negative binomial with this mean and size (variance mean + mean^2 / size);
    Poisson if size is inf."""
    h = np.asarray(h, dtype=np.int64)
    mean, size = np.broadcast_to(mean, h.shape).astype(float), np.broadcast_to(size, h.shape).astype(float)
    out = np.ones(h.shape)
    k = h > 0
    poi = k & np.isinf(size)
    out[poi] = gammainc(h[poi], mean[poi])
    nb = k & ~np.isinf(size)
    out[nb] = betainc(h[nb], size[nb], mean[nb] / (size[nb] + mean[nb]))
    return out


def profile(pos, length, tx, anchor, lengths, window=PROFILE):
    """P-sites per (read length in [lo, hi], nt from the anchor of their transcript in window). anchor: a flat
    position per transcript (-1: none), e.g. its annotated start."""
    lo, hi = lengths
    t = tx[pos]
    c = anchor[np.maximum(t, 0)]
    d = pos - c
    k = (t >= 0) & (c >= 0) & (d >= window[0]) & (d < window[1]) & (length >= lo) & (length <= hi)
    n = window[1] - window[0]
    return np.bincount((length[k] - lo) * n + d[k] - window[0], minlength=(hi - lo + 1) * n).reshape(-1, n)


def region_counts(pos, length, region, lengths):
    """P-sites per (read length in [lo, hi], region: 5' UTR, CDS, 3' UTR, transcript)."""
    lo, hi = lengths
    r = region[pos]
    k = (r >= 0) & (length >= lo) & (length <= hi)
    n_reg = len(REGIONS)
    return np.bincount((length[k] - lo) * n_reg + r[k], minlength=(hi - lo + 1) * n_reg).reshape(-1, n_reg)


def cds_bins(pos, length, tx, layout, lengths):
    """P-sites per (read length in [lo, hi], bin of CDS_BINS along the annotated CDS, frame of the CDS)."""
    lo, hi = lengths
    t = np.maximum(tx[pos], 0)
    cs, ce = layout.cds_start[t], layout.cds_end[t]
    x = pos - layout.first[t] - cs
    k = (tx[pos] >= 0) & (cs >= 0) & (x >= 0) & (x < ce - cs) & (length >= lo) & (length <= hi)
    b = x[k] * CDS_BINS // (ce - cs)[k]
    return np.bincount(((length[k] - lo) * CDS_BINS + b) * 3 + x[k] % 3,
                       minlength=(hi - lo + 1) * CDS_BINS * 3).reshape(-1, CDS_BINS, 3)


def interior_frames(pos, length, tx, layout, lengths):
    """P-sites per (read length in [lo, hi], frame of the annotated CDS) at least INTERIOR nt inside the CDS at both
    ends: no start or stop peaks."""
    lo, hi = lengths
    t = np.maximum(tx[pos], 0)
    cs, ce = layout.cds_start[t], layout.cds_end[t]
    x = pos - layout.first[t]
    k = (tx[pos] >= 0) & (cs >= 0) & (x >= cs + INTERIOR) & (x < ce - INTERIOR) & (length >= lo) & (length <= hi)
    return np.bincount((length[k] - lo) * 3 + (x[k] - cs[k]) % 3, minlength=(hi - lo + 1) * 3).reshape(-1, 3)


class LengthCounts(NamedTuple):
    """One library's P-sites per read length."""
    start: np.ndarray     # at nt PROFILE from the annotated starts (profile())
    frames: np.ndarray    # per CDS frame, in the CDS interior (interior_frames())
    regions: np.ndarray   # per region (region_counts())


def peak_beside(prof, allowed):
    """Per row of prof: the column of the peak among the allowed columns (allowed: (rows, columns) bool; the first if
    tied), its value (0 if no column is allowed), and the values in the columns before and after it (0 off the
    profile)."""
    n = np.arange(len(prof))
    i = np.argmax(np.where(allowed, prof, -1), axis=1)
    at = np.where(allowed[n, i], prof[n, i], 0)
    padded = np.pad(prof, ((0, 0), (1, 1)))
    return i, at, padded[n, i], padded[n, i + 2]


def calibrate(prof, frames, lengths, pair=None):
    """Per read length of one library (prof: P-sites at nt PROFILE from the annotated starts; frames: interior P-sites
    per CDS frame): the start peak, the frames, whether the length is calibrated, and its shift (the P-sites move by
    it). pair: per read length, True for a pair length (default: none). The anchor is the start peak, or for a pair
    length the pair of the peak and the larger of the two nt beside it (beside_nt: the nt before the peak if tied), at
    its first nt. Calibrated: the peak and the anchor are within MAX_SHIFT nt of nt 0, the peak has MIN_ANCHOR P-sites
    or more, and the frame is not clearly against the anchor (the dominant frame, with MIN_FRAME interior P-sites or
    more, is not a frame of the anchor and has FRAME_GUARD times the P-sites of the larger frame that the shift moves
    onto frame 0 or, for a pair, 1). Shift: minus the anchor; for a length that is not calibrated, the shift that
    moves a clear dominant frame onto frame 0 (frame 1: -1, frame 2: +1), or for a pair length a clear frame pair (the
    two consecutive frames with the most P-sites, the smaller with FRAME_GUARD times the third frame) onto frames 0 and
    1, else 0. Mixed: MIN_ANCHOR P-sites or more at the peak, and the larger of the two nt beside it (beside_peak, at
    nt beside_nt) has MIXED of them or more (beside_share: beside_peak / at_peak)."""
    n = np.arange(len(prof))
    pair = np.zeros(len(prof), dtype=bool) if pair is None else np.asarray(pair, dtype=bool)
    to_kernel_end = np.arange(prof.shape[1]) < KERNEL_END - PROFILE[0]
    col, at_peak, before, after = peak_beside(prof, np.broadcast_to(to_kernel_end, prof.shape))
    peak = col + PROFILE[0]
    beside = np.maximum(before, after)
    beside_nt = np.where(before >= after, peak - 1, peak + 1)
    mixed = (at_peak >= MIN_ANCHOR) & (beside >= MIXED * at_peak)
    anchor = np.where(pair, np.minimum(peak, beside_nt), peak)
    total = frames.sum(axis=1)
    dominant = frames.argmax(axis=1)
    own = (dominant == anchor % 3) | (pair & (dominant == (anchor + 1) % 3))
    held = np.where(pair, np.maximum(frames[n, anchor % 3], frames[n, (anchor + 1) % 3]), frames[n, anchor % 3])
    against = (total >= MIN_FRAME) & ~own & (frames[n, dominant] >= FRAME_GUARD * held)
    # the peak too: a codon-1 peak at +3 with its neighbour at +2 is not a pair at the start codon
    calibrated = (np.abs(peak) <= MAX_SHIFT) & (np.abs(anchor) <= MAX_SHIFT) & (at_peak >= MIN_ANCHOR) & ~against
    top2 = np.sort(frames, axis=1)[:, -2:]
    clear = (total >= MIN_FRAME) & (top2[:, 1] >= FRAME_GUARD * top2[:, 0])
    duo = (frames + np.roll(frames, -1, axis=1)).argmax(axis=1)   # the frame pair (duo, duo + 1) with the most
    clear_duo = (total >= MIN_FRAME) & (np.minimum(frames[n, duo], frames[n, (duo + 1) % 3])
                                        >= FRAME_GUARD * frames[n, (duo + 2) % 3])
    by_frame = np.where(pair, np.where(clear_duo, np.array([0, -1, 1])[duo], 0),
                        np.where(clear, np.array([0, -1, 1])[dominant], 0))
    shift = np.where(calibrated, -anchor, by_frame)
    return pd.DataFrame({"length": np.arange(lengths[0], lengths[1] + 1), "start_peak": peak, "at_peak": at_peak,
                         "beside_peak": beside, "beside_nt": beside_nt, "beside_share": beside / np.maximum(at_peak, 1),
                         "mixed": mixed, "pair": pair,
                         "frame0": frames[:, 0], "frame1": frames[:, 1], "frame2": frames[:, 2],
                         "calibrated": calibrated, "shift": shift})


def pair_lengths(mixed):
    """Per read length (rows: the harringtonine libraries' calibrate() rows, columns length and mixed): a pair length
    is mixed in half of these libraries or more."""
    m = mixed.groupby("length")["mixed"]
    return (2 * m.sum() >= m.size()).sort_index().to_numpy()


def length_qc(prof, frames, regions, lengths, pair=None):
    """Per read length of one library: calibrate()'s columns, the P-sites (all, their share of the library's, and per
    region), the top frame's share of the CDS-interior P-sites, and whether the length is periodic (MIN_FRAME
    interior P-sites or more, MIN_TOP of them in the top frame)."""
    out = calibrate(prof, frames, lengths, pair)
    psites = regions.sum(axis=1)
    interior = frames.sum(axis=1)
    top = frames.max(axis=1) / np.maximum(interior, 1)
    out["psites"] = psites
    out["share"] = psites / max(psites.sum(), 1)
    for i, r in enumerate(REGIONS):
        out[r] = regions[:, i]
    out["top_frame"] = top
    out["periodic"] = (interior >= MIN_FRAME) & (top >= MIN_TOP)
    return out


def choose_start_lengths(table, window):
    """The QC verdict per read length ("pass", or why not) and the start lengths. table: length_qc's rows of every
    library, with kind. A start length is in the window, calibrated in every harringtonine library and periodic in
    every library; of the runs of consecutive such lengths, the one with the most P-sites (the first if tied)."""
    harr = table["kind"] == "harringtonine"
    n_harr, n_lib = table.loc[harr, "library"].nunique(), table["library"].nunique()
    uncal = (~table.loc[harr, "calibrated"]).groupby(table.loc[harr, "length"]).sum()
    aperiodic = (~table["periodic"]).groupby(table["length"]).sum()
    reads = table.groupby("length")["psites"].sum()
    why = {}
    for length in reads.index:
        r = []
        if not window[0] <= length <= window[1]:
            r.append(f"outside the window {window[0]}-{window[1]}")
        if uncal.get(length, 0):
            r.append(f"not calibrated in {uncal[length]} of {n_harr} harringtonine libraries")
        if aperiodic[length]:
            r.append(f"no clear frame in {aperiodic[length]} of {n_lib} libraries")
        why[int(length)] = "; ".join(r)
    runs = []
    for length in (x for x in why if not why[x]):
        if runs and runs[-1][-1] == length - 1:
            runs[-1].append(length)
        else:
            runs.append([length])
    best = max(runs, key=lambda run: reads[run].sum()) if runs else []
    for length in (x for x in why if not why[x]):
        why[length] = "pass" if length in best else f"outside the run {best[0]}-{best[-1]}"
    return pd.Series(why), best


def qc_profiles(start, stop, cds, lo):
    """One library's profiles as a long table: P-sites per read length at nt x from the annotated starts (anchor
    start) and from the first nt of the stop codons (stop), and per bin x along the CDSs (cds); frame: of the CDS."""
    rows = []
    for anchor, arr, x0 in (("start", start, PROFILE[0]), ("stop", stop, STOP_PROFILE[0])):
        n_len, n = arr.shape
        x = np.tile(np.arange(x0, x0 + n), n_len)
        rows.append(pd.DataFrame({"length": np.repeat(np.arange(lo, lo + n_len), n), "anchor": anchor, "x": x,
                                  "frame": x % 3, "psites": arr.ravel()}))
    n_len = cds.shape[0]
    rows.append(pd.DataFrame({"length": np.repeat(np.arange(lo, lo + n_len), CDS_BINS * 3), "anchor": "cds",
                              "x": np.tile(np.repeat(np.arange(CDS_BINS), 3), n_len),
                              "frame": np.tile(np.arange(3), n_len * CDS_BINS), "psites": cds.ravel()}))
    return pd.concat(rows, ignore_index=True)


def shift_psites(psites, shift, lo, tx):
    """P-sites (PSites) moved by the shift of their length (shift[length - lo]); those that leave their transcript are
    dropped."""
    pos, length, multi = psites
    new = pos + shift[length - lo]
    k = (new >= 0) & (new < len(tx))
    k[k] = tx[new[k]] == tx[pos[k]]
    return PSites(new[k], length[k], multi[k])


def kernel(hp, ep):
    """Per read length: harringtonine minus elongation P-sites (scaled on codons 11-20), as a share of
    nt -15..32, from profiles at the annotated starts."""
    far = slice(KERNEL_END - PROFILE[0], PROFILE[1] - PROFILE[0])
    k = hp[:, far].sum(axis=1) / np.maximum(ep[:, far].sum(axis=1), 1)
    kern = np.maximum(hp - k[:, None] * ep, 0)[:, :KERNEL_END - PROFILE[0]]
    return kern / np.maximum(kern.sum(axis=1, keepdims=True), 1e-300)


@dataclass
class Pool:
    h: np.ndarray           # harringtonine P-sites of the start lengths per position, summed over the pool
    e: np.ndarray           # elongation
    multi: np.ndarray       # multimapped P-sites, both kinds
    f: np.ndarray = None    # factor per bin
    s: np.ndarray = None    # scale per transcript
    dens: np.ndarray = None  # prior mean elongation reads per nt, per (transcript, zone)
    a: np.ndarray = None    # prior shape per region (inf: no spread)


def moments_nb(h, m, base_var):
    """Moment estimate of phi in var = base_var + phi m^2, at least 0."""
    return float(max((((h - m) ** 2).sum() - base_var.sum()) / max((m ** 2).sum(), 1e-300), 0))


def fit_pool(pool, tx, region, zone, bins, bg, null_pos, null_win, n_tx):
    """The pool's background: factors and scales (fit_factors), the prior mean and shape of the
    elongation reads of a window."""
    bg_pos = np.flatnonzero(bg)
    pool.f, pool.s = fit_factors(pool.h, pool.e, tx, bins, null_pos, bg_pos, n_tx, len(bin_names()), PRIOR_READS)
    key = tx[bg_pos].astype(np.int64) * ZONES + zone[bg_pos]
    z_mean = (np.bincount(zone[bg_pos], weights=pool.e[bg_pos], minlength=ZONES)
              / np.maximum(np.bincount(zone[bg_pos], minlength=ZONES), 1))
    pool.dens = ((np.bincount(key, weights=pool.e[bg_pos], minlength=n_tx * ZONES) + PRIOR_NT * np.tile(z_mean, n_tx))
                 / (np.bincount(key, minlength=n_tx * ZONES) + PRIOR_NT))
    e = window_sums(pool.e, null_win)
    lam0 = prior_mean(pool, null_win, tx, zone)
    pool.a = np.full(len(REGIONS), np.inf)
    for r in range(len(REGIONS)):
        k = region[null_win] == r
        inv = moments_nb(e[k].astype(float), lam0[k], lam0[k]) if k.any() else 0.0
        pool.a[r] = 1 / inv if inv > 0 else np.inf


def prior_mean(pool, pos, tx, zone):
    """Prior mean elongation reads of the windows at pos: the densities of their nt's zones in their transcript."""
    t = tx[pos].astype(np.int64)
    lam = np.zeros(len(pos))
    for k in range(WIN[0], WIN[1]):
        z = zone[pos + k]
        lam += np.where((z >= 0) & (tx[pos + k] == t), pool.dens[t * ZONES + np.maximum(z, 0)], 0.0)
    return lam


class Moments(NamedTuple):
    """Per window: its harringtonine and elongation reads, and the mean and variance (without phi) of its
    harringtonine reads."""
    h: np.ndarray
    e: np.ndarray
    mean: np.ndarray
    var: np.ndarray


def window_moments(pool, pos, tx, region, zone, bins):
    """Harringtonine reads in the windows at pos, and their mean and variance without phi: the
    elongation rate of a window has a gamma prior (mean lam0, shape a), updated by its e reads. The
    mean is R times the larger of e and the posterior mean: the prior lifts a window with few elongation
    reads, but does not shrink an elongation peak (a pause)."""
    t, r = tx[pos].astype(np.int64), region[pos]
    e = window_sums(pool.e, pos)
    R = pool.f[bins[pos]] * pool.s[t]
    lam0 = prior_mean(pool, pos, tx, zone)
    a = pool.a[r]
    fin = np.isfinite(a)
    af = np.where(fin, a, 0.0)
    m = R * np.maximum(e, np.where(fin, (af + e) / (af / lam0 + 1), lam0))
    var = m + np.where(fin, m ** 2 / np.where(fin, af + e, 1), 0.0)
    return Moments(window_sums(pool.h, pos), e, m, var)


def pooled_moments(pools, pos, tx, region, zone, bins):
    """window_moments() of the windows at pos, summed over the pools: pools add their means and variances."""
    got = [window_moments(pl, pos, tx, region, zone, bins) for pl in pools.values()]
    return Moments(*(sum(col) for col in zip(*got)))


def nb_size(m, var, phi):
    """NB size with mean m and variance var + phi m^2; inf (Poisson) if that is not above m."""
    extra = var - m + phi * m ** 2
    return np.where(extra > 0, m ** 2 / np.maximum(extra, 1e-300), np.inf)


def tail_phi(h, m, var):
    """The smallest phi on PHI_GRID with at most TAIL of the windows at p <= TAIL."""
    for phi in PHI_GRID:
        if (nb_sf(h, m, nb_size(m, var, phi)) <= TAIL).mean() <= TAIL:
            return float(phi)
    return float(PHI_GRID[-1])


def depth_bin(m):
    return np.searchsorted(DEPTH, m, "right") - 1


def calibrate_phi(h, m, var, region):
    """phi per (region, depth bin) from null windows; MIN_NULL or more windows each, else the region's,
    else all windows'."""
    d = depth_bin(m)
    every = tail_phi(h, m, var)
    phi = {}
    for r in range(len(REGIONS)):
        k = region == r
        own = tail_phi(h[k], m[k], var[k]) if k.sum() >= MIN_NULL else every
        for b in range(len(DEPTH)):
            kb = k & (d == b)
            phi[r, b] = tail_phi(h[kb], m[kb], var[kb]) if kb.sum() >= MIN_NULL else own
    return phi


def phi_stats(phi):
    """phi as starts_stats.tsv rows."""
    return {f"phi {REGIONS[r]} {DEPTH[b]}": value for (r, b), value in phi.items()}


def window_phi(m, region, phi):
    return np.array([phi[r, b] for r, b in zip(region.tolist(), depth_bin(m).tolist())], dtype=float)


def pvalues(h, m, var, region, phi):
    """p of each window: h ~ NB(m, var + phi m^2) with the phi of its region and depth bin."""
    return nb_sf(h, m, nb_size(m, var, window_phi(m, region, phi)))


def given_starts(path, layout):
    """The starts of an ORF table (TSV ORF_id Name start), plus the annotated starts: flat position ->
    ORF_ids. Rows whose transcript is not in the layout or whose start codon is off it are counted."""
    table = pd.read_csv(path, sep="\t", dtype={"ORF_id": str, "Name": str})[["ORF_id", "Name", "start"]]
    idx = pd.Series(np.arange(len(layout.names)), index=layout.names)
    t = table["Name"].map(idx)
    t, s = t.fillna(0).to_numpy(np.int64), table["start"].to_numpy(np.int64)
    ok = table["Name"].isin(idx.index).to_numpy() & (s >= 0) & (s + 3 <= layout.tx_len[t])
    c = np.flatnonzero(layout.cds_start >= 0)
    ids = pd.concat([pd.Series(table["ORF_id"].to_numpy()[ok], index=layout.first[t[ok]] + s[ok]),
                     pd.Series([f"{n}:CDS" for n in layout.names[c]], index=layout.first[c] + layout.cds_start[c])])
    ids = ids.groupby(level=0, sort=True).agg(",".join)
    return ids.index.to_numpy(np.int64), ids.to_numpy(), {"orf_table_rows": len(table),
                                                          "orf_table_rows_dropped": int((~ok).sum())}


def scan_candidates(layout, tx, pools, cs_pos, codons=START_CODONS):
    """Flat positions of the start codons (`codons`) with harringtonine reads in their window, and the
    annotated starts."""
    c = np.flatnonzero(codon_mask(layout, codons) & (tx >= 0))
    h = sum(window_sums(pl.h, c) for pl in pools.values())
    return np.union1d(c[h > 0], cs_pos)


def next_stops(layout):
    """Per flat position x: the first in-frame stop codon at x or after (same frame as x), or the array size."""
    stop = codon_mask(layout, tuple(sorted(annotation.STOPS)))
    out = np.empty(layout.size, dtype=np.int64)
    for r in range(3):
        i = np.arange(r, layout.size, 3)
        out[i] = np.minimum.accumulate(np.where(stop[i], i, layout.size)[::-1])[::-1]
    return out


def orf_types(start, end, cs, ce):
    """ORF type from transcript coordinates against the annotated CDS [cs, ce) (cs -1: no CDS)."""
    up, inside = start < cs, (start > cs) & (start < ce)
    frame0 = (start - cs) % 3 == 0
    return np.select([cs < 0, start == cs, up & (end <= cs), up & frame0, up, inside & frame0, inside],
                     ["other", "CDS", "uORF", "extension", "uoORF", "truncation", "internal"], "dORF")


def bh(p, m):
    """Benjamini-Hochberg q of p among m tests; the m - len(p) tests not given have p = 1."""
    order = np.argsort(p, kind="stable")
    q = np.minimum.accumulate((p[order] * m / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(q, 1.0)
    return out


def near_starts(pos, p, passing, annotated=()):
    """Called starts. The candidates in codons 1-10 of an annotated start (indices `annotated`) are not called:
    harringtonine raises them when it is used. Then by p (then position), a passing candidate that no called start
    stops is called. It stops the candidates within NEAR nt and in its codons 1-10. Returns called, and per
    candidate the index of the called or annotated start that stops it (-1: none)."""
    called = np.zeros(len(pos), dtype=bool)
    stopped_by = np.full(len(pos), -1, dtype=np.int64)
    for i in annotated:
        lo, hi = np.searchsorted(pos, [pos[i] + 1, pos[i] + RAISED[1]])
        stopped_by[lo:hi] = i
    for i in np.lexsort((pos, p)):
        if not passing[i] or stopped_by[i] >= 0:
            continue
        called[i] = True
        lo, hi = np.searchsorted(pos, [pos[i] + NEAR[0], pos[i] + RAISED[1]])
        near = np.arange(lo, hi)
        near = near[~called[near] & (stopped_by[near] < 0)]
        stopped_by[near] = i
    return called, stopped_by


def codon_table(layout, tested, pos, cls, called):
    """Per class (CLASSES) and start codon: the tests (tested: flat positions per class) and the calls (pos[called],
    their classes cls[called]), each also as a share of the class's. Codons in the order of START_CODONS, then the
    others."""
    def counts(x):
        s = layout.seq
        key = s[x].astype(np.int64) * 65536 + s[x + 1].astype(np.int64) * 256 + s[x + 2]
        u, n = np.unique(key, return_counts=True)
        return {bytes([k >> 16, (k >> 8) & 255, k & 255]).decode(): int(c) for k, c in zip(u.tolist(), n.tolist())}
    rows = []
    for c, name in enumerate(CLASSES):
        nt, nc = counts(tested[c]), counts(pos[called & (cls == c)])
        st, sc = sum(nt.values()), sum(nc.values())
        order = [x for x in START_CODONS if x in nt or x in nc]
        order += sorted((set(nt) | set(nc)) - set(order))
        rows += [[name, x, nt.get(x, 0), nt.get(x, 0) / max(st, 1), nc.get(x, 0), nc.get(x, 0) / max(sc, 1)]
                 for x in order]
    return pd.DataFrame(rows, columns=CODON_COLUMNS)


def class_tests(layout, tx, region, pos, cs_pos, scan, codons):
    """The BH class (CLASSES) of each candidate at pos, and per class the number of its tests and their flat positions.
    A scan tests every start codon (of `codons`) of the class, also those without harringtonine reads (p = 1)."""
    cls = np.where(np.isin(pos, cs_pos), ANNOTATED, region[pos])
    if scan:
        every = np.flatnonzero(codon_mask(layout, codons) & (tx >= 0))
        every = every[~np.isin(every, cs_pos)]
        n_tests = np.r_[np.bincount(region[every], minlength=len(REGIONS)), len(cs_pos)]
        tested = [every[region[every] == r] for r in range(len(REGIONS))] + [cs_pos]
    else:
        n_tests = np.bincount(cls, minlength=len(CLASSES))
        tested = [pos[cls == c] for c in range(len(CLASSES))]
    return cls, n_tests, tested


def class_q(p, cls, n_tests):
    """BH q per class: the candidates of class c (cls == c) among the n_tests[c] tests of the class."""
    q = np.full(len(p), np.nan)
    for c in range(len(CLASSES)):
        k = cls == c
        if k.any():
            q[k] = bh(p[k], n_tests[c])
    return q


def class_stats(p, q, cls, n_tests, null_p, null_region, fdr):
    """Per class: its tests, candidates and calls (q < fdr); the share of the null windows at or below the p of its
    weakest call (the null windows of the class's region; for the annotated starts, the CDS's), and the false calls
    that share expects in its tests."""
    stats = {}
    for c, name in enumerate(CLASSES):
        k = cls == c
        if not k.any():
            continue
        calls = q[k] < fdr
        cut = p[k][calls].max() if calls.any() else 0.0
        kn = null_region == (CDS if c == ANNOTATED else c)
        rate = (null_p[kn] <= cut).mean() if kn.any() and cut > 0 else 0.0
        stats.update({f"tests {name}": int(n_tests[c]), f"candidates {name}": int(k.sum()),
                      f"q<{fdr:g} {name}": int(calls.sum()), f"null_rate_at_cut {name}": rate,
                      f"expected_false {name}": rate * n_tests[c]})
    return stats


def scan_calls(out, pos, p, q, cls, layout, tx, fdr, stats):
    """A scan's calls (near_starts() on q < fdr) and the ORF of each candidate to its first in-frame stop, typed
    against the annotated CDS: adds called, stopped_by, end and type to out (the starts table), and the counts to
    stats. Returns the ORFs of the called starts but the annotated ones (start_orfs.tsv), and the catalogue
    (catalogue.tsv): the uORFs and uoORFs of the called starts, and every annotated CDS (end: the annotation's)."""
    t = tx[pos].astype(np.int64)
    start = pos - layout.first[t]
    called, stopped_by = near_starts(pos, p, q < fdr, np.flatnonzero(cls == ANNOTATED))
    stop = next_stops(layout)[pos + 3]
    has_stop = (stop < layout.size) & (tx[np.minimum(stop, layout.size - 1)] == tx[pos])
    end = np.where(has_stop, stop - layout.first[t], -1)
    out["called"] = called
    out["stopped_by"] = pd.Series(out["ORF_id"].to_numpy()[np.maximum(stopped_by, 0)]).where(stopped_by >= 0)
    out["end"] = pd.Series(end, dtype="Int64").mask(~has_stop)
    out["type"] = np.where(has_stop, orf_types(start, end, layout.cds_start[t], layout.cds_end[t]), "no stop")
    typ = out["type"].to_numpy()
    orfs = out[called & has_stop & (typ != "CDS")].rename(columns={"codon": "start_codon"})[ORF_COLUMNS]
    cat = out[(called & np.isin(typ, ["uORF", "uoORF"])) | (typ == "CDS")].rename(columns={"codon": "start_codon"})
    cat = cat.assign(end=np.where(cat["type"] == "CDS", layout.cds_end[t[cat.index]], cat["end"]))[ORF_COLUMNS]
    stats.update({f"catalogue {k}": int(v) for k, v in cat["type"].value_counts().sort_index().items()})
    stats.update(called=int(called.sum()), called_without_stop=int((called & ~has_stop).sum()),
                 stopped=int((stopped_by >= 0).sum()))
    stats.update({f"start_orfs {k}": int(v) for k, v in orfs["type"].value_counts().sort_index().items()})
    return orfs, cat


def shared_transcripts(prefixes):
    """The transcripts in the orfs.tsv of every orfs run."""
    names = None
    for prefix in prefixes:
        n = set(pd.read_csv(f"{prefix}.orfs.tsv", sep="\t", usecols=["Name"], dtype=str)["Name"])
        names = n if names is None else names & n
    return names


def length_tables(psites, role, pool_of, layout, tx, region, cs_flat, ce_flat, lengths):
    """The QC of each read length in each library (lengths.tsv but its columns qc and start_length) and the profiles
    behind it (profiles.tsv), from the P-sites at the pool's offsets. psites, role, pool_of: per library."""
    counts, profiles = {}, []
    for prefix, (pos, length, _) in psites.items():
        start = profile(pos, length, tx, cs_flat, lengths)
        counts[prefix] = LengthCounts(start, interior_frames(pos, length, tx, layout, lengths),
                                      region_counts(pos, length, region, lengths))
        profiles.append(qc_profiles(start, profile(pos, length, tx, ce_flat, lengths, STOP_PROFILE),
                                    cds_bins(pos, length, tx, layout, lengths), lengths[0]).assign(library=prefix))
    # one anchor per read length for all libraries: a pair if the length is mixed in half of the harringtonine
    # libraries or more
    pair = pair_lengths(pd.concat([calibrate(c.start, c.frames, lengths) for prefix, c in counts.items()
                                   if role[prefix] == "harringtonine"]))
    qc = []
    for prefix, c in counts.items():
        q = length_qc(c.start, c.frames, c.regions, lengths, pair).assign(library=prefix, pool=pool_of[prefix],
                                                                         kind=role[prefix])
        orfs_offset = pd.read_csv(f"{prefix}.frames.tsv", sep="\t").set_index("length")["offset"].astype("Int64")
        q["offset"] = q["length"].map(orfs_offset).astype("Int64") + q["shift"]   # the offset the P-sites end at
        qc.append(q)
    table = pd.concat(qc, ignore_index=True)
    interior = table[["frame0", "frame1", "frame2"]].sum(axis=1)
    table["interior_share"] = interior / interior.groupby(table["library"]).transform("sum").clip(lower=1)
    return table, pd.concat(profiles, ignore_index=True)[["library", "length", "anchor", "x", "frame", "psites"]]


def pair_distances(sheet, table):
    """pairs.tsv: per pair (a row of the sheet), the read-length distance: half the summed absolute difference of the
    two libraries' interior_share (table: length_tables())."""
    interior = table.pivot(index="length", columns="library", values="interior_share").fillna(0)
    pairs = sheet[["pool", "harringtonine", "elongation"]].copy()
    pairs["distance"] = [0.5 * float((interior[h] - interior[e]).abs().sum())
                         for h, e in zip(pairs["harringtonine"], pairs["elongation"])]
    return pairs[PAIR_COLUMNS]


def library_table(psites, start_lengths, table, layout, tx, region, cs_pos, role, pool_of, lengths):
    """libraries.tsv: per library, its start-length P-sites after the calibration, by region and by frame inside the
    CDSs; the start peak, per start, over the P-sites per in-frame codon inside the CDSs (the codons that
    interior_frames counts in frame 0). Not the stop peak: on real data it did not fade with harringtonine (7.5-20x
    against 8-18x)."""
    has = layout.cds_start >= 0
    codons = int(((np.maximum(layout.cds_end[has] - layout.cds_start[has] - 2 * INTERIOR, 0) + 2) // 3).sum())
    rows = []
    for prefix, (pos, length, _) in psites.items():
        k = np.isin(length, start_lengths)
        reg = np.bincount(region[pos[k]], minlength=len(REGIONS))
        fr = interior_frames(pos[k], length[k], tx, layout, lengths).sum(axis=0)
        total = int(table.loc[table["library"] == prefix, "psites"].sum())
        at_start = int(np.isin(pos[k], cs_pos).sum())
        per_codon = fr[0] / max(codons, 1)
        rows.append({"library": prefix, "pool": pool_of[prefix], "kind": role[prefix], "psites": total,
                     "start_psites": int(k.sum()), "start_share": k.sum() / max(total, 1),
                     **dict(zip(REGIONS, reg / max(reg.sum(), 1))), "interior": int(fr.sum()),
                     **{f"frame{f}": x for f, x in enumerate(fr / max(fr.sum(), 1))}, "at_start": at_start,
                     "start_enrichment": (at_start / len(cs_pos) / per_codon if per_codon > 0 and len(cs_pos)
                                          else np.nan)})
    return pd.DataFrame(rows, columns=LIBRARY_COLUMNS)


def kernel_table(sheet, psites, tx, cs_flat, lengths):
    """kernel.tsv: per pool and read length, the harringtonine and elongation P-sites (after the calibration) at nt
    -15..32 from the annotated starts, and the kernel."""
    lo, hi = lengths
    prof = {prefix: profile(pos, length, tx, cs_flat, lengths) for prefix, (pos, length, _) in psites.items()}
    nt = np.arange(PROFILE[0], KERNEL_END)
    rows = []
    for name, group in sheet.groupby("pool", sort=True):
        hp = sum(prof[prefix] for prefix in group["harringtonine"])
        ep = sum(prof[prefix] for prefix in group["elongation"])
        for length, kr in zip(range(lo, hi + 1), kernel(hp, ep)):
            i = length - lo
            rows.append(pd.DataFrame({"pool": name, "length": length, "nt": nt,
                                      "harringtonine": hp[i, :len(nt)], "elongation": ep[i, :len(nt)],
                                      "kernel": kr}))
    return pd.concat(rows, ignore_index=True)


def make_pools(sheet, psites, start_lengths, size):
    """Per pool (by name): the start-length P-sites per flat position, summed over its libraries (Pool)."""
    pools = {}
    for name, group in sheet.groupby("pool", sort=True):
        sums = {kind: np.zeros(size, dtype=np.int64) for kind in ("harringtonine", "elongation", "multi")}
        for kind in ("harringtonine", "elongation"):
            for prefix in group[kind]:
                pos, length, multi = psites[prefix]
                k = np.isin(length, start_lengths)
                sums[kind] += np.bincount(pos[k], minlength=size)
                sums["multi"] += np.bincount(pos[k & multi], minlength=size)
        pools[name] = Pool(sums["harringtonine"], sums["elongation"], sums["multi"])
    return pools


def null_windows(pools, nulls, layout, tx):
    """The null codons with mostly unique reads (multimapped share below MAX_MULTI), and of those the null windows:
    one null codon per STRIDE nt of a transcript, so the windows do not overlap."""
    reads = sum(window_sums(pl.h, nulls) + window_sums(pl.e, nulls) for pl in pools.values())
    multi = sum(window_sums(pl.multi, nulls) for pl in pools.values())
    null_pos = nulls[multi < MAX_MULTI * np.maximum(reads, 1)]
    null_win = null_pos[(null_pos - layout.first[tx[null_pos]]) % STRIDE == 0]
    return null_pos, null_win


def finish(stats, out_prefix, report):
    """starts_stats.tsv, then the report (report/starts.py) if `report`."""
    write_stats(stats, f"{out_prefix}.starts_stats.tsv")
    if report:
        from . import report as rpt
        rpt.starts_report(out_prefix)


def run(sheet_path, gtf, fasta, out_prefix, orfs_path=None, scan=False, start_lengths=None, start_codons=None,
        fdr=FDR, window=WINDOW, report=True):
    """Start evidence from the `pitstop orfs` runs of a library sheet (TSV pool, harringtonine, elongation:
    out prefixes, one row per replicate), for the starts of an ORF table and the annotated starts, or for every
    start codon (of `start_codons`, default START_CODONS) with harringtonine reads (scan: also the starts called
    at q < fdr and their ORFs). start_lengths (lo, hi): the read lengths that count; default: the QC's, in `window`.
    report: also <out_prefix>.starts.html (report/starts.py), from the tables only; also when no read length passes the QC,
    before the error."""
    if (orfs_path is None) == (not scan):
        raise ValueError("give an ORF table or scan, not both")
    if start_codons is not None and not scan:
        raise ValueError("start codons choose the scan's candidates: give them with a scan")
    start_codons = START_CODONS if start_codons is None else tuple(start_codons)
    # the inputs, for the report: file names only, no paths
    stats = {"pitstop": __version__, "libraries_sheet": os.path.basename(sheet_path), "gtf": os.path.basename(gtf),
             "fasta": os.path.basename(fasta), "candidates_from": "scan" if scan else os.path.basename(orfs_path),
             "fdr": fdr, "length_window": f"{window[0]}-{window[1]}",
             "start_lengths_given": "" if start_lengths is None else f"{start_lengths[0]}-{start_lengths[1]}"}
    sheet = pd.read_csv(sheet_path, sep="\t", dtype=str)
    prefixes = list(sheet["harringtonine"]) + list(sheet["elongation"])
    role, pool_of = {}, {}   # per library: its kind and its pool
    for kind in ("harringtonine", "elongation"):
        for prefix, pool in zip(sheet[kind], sheet["pool"]):
            role.setdefault(prefix, kind)
            pool_of.setdefault(prefix, pool)

    layout = make_layout(gtf, fasta, shared_transcripts(prefixes))
    tx, region, dist = layout.positions()
    bins = bin_index(region, dist)
    zone = zone_index(layout, tx, region)
    n_tx = len(layout.names)
    bg = background(layout, tx)
    nulls = null_codons(bg, codon_mask(layout))
    cs_pos, cs_tx = layout.cds_starts()
    cs_flat, ce_flat = layout.flat(layout.cds_start), layout.flat(layout.cds_end)
    stats.update(pools=sheet["pool"].nunique(), libraries=len(prefixes), transcripts=n_tx,
                 transcripts_with_cds=len(cs_tx), null_codons=len(nulls))

    # the QC of the read lengths, from the P-sites at the pool's offsets: the start lengths and the shifts
    psites = {prefix: read_psites(prefix, layout, (0, np.iinfo(np.int16).max)) for prefix in prefixes}
    lo = min(int(x.length.min()) for x in psites.values())
    hi = max(int(x.length.max()) for x in psites.values())
    table, profiles = length_tables(psites, role, pool_of, layout, tx, region, cs_flat, ce_flat, (lo, hi))
    verdict, chosen = choose_start_lengths(table, window)
    start_lengths = chosen if start_lengths is None else list(range(start_lengths[0], start_lengths[1] + 1))
    table["qc"] = table["length"].map(verdict)
    table["start_length"] = table["length"].isin(start_lengths)
    write_tsv(table[LENGTH_COLUMNS], f"{out_prefix}.lengths.tsv")
    write_tsv(pair_distances(sheet, table), f"{out_prefix}.pairs.tsv")
    write_tsv(profiles, f"{out_prefix}.profiles.tsv")
    del profiles
    if not start_lengths:
        stats["start_lengths"] = ""
        finish(stats, out_prefix, report)
        raise ValueError(f"no read length passes the QC (in the window {window[0]}-{window[1]}, calibrated in "
                         f"every harringtonine library, a clear frame in every library): see "
                         f"{out_prefix}.lengths.tsv, or choose them with --start-lengths")
    stats["start_lengths"] = ",".join(map(str, start_lengths))
    log.info("start lengths: %s", start_lengths)

    # the calibrated P-sites, per library and summed per pool
    shift = {prefix: g["shift"].to_numpy() for prefix, g in table.groupby("library", sort=False)}
    psites = {prefix: shift_psites(x, shift[prefix], lo, tx) for prefix, x in psites.items()}
    write_tsv(library_table(psites, start_lengths, table, layout, tx, region, cs_pos, role, pool_of, (lo, hi)),
              f"{out_prefix}.libraries.tsv")
    kernels = kernel_table(sheet, psites, tx, cs_flat, (lo, hi))
    pools = make_pools(sheet, psites, start_lengths, layout.size)
    del psites
    for name, pl in pools.items():
        stats[f"harringtonine_psites {name}"] = int(pl.h.sum())
        stats[f"elongation_psites {name}"] = int(pl.e.sum())

    # the background of each pool and phi, from the null windows
    null_pos, null_win = null_windows(pools, nulls, layout, tx)
    stats.update(null_codons_unique=len(null_pos), null_windows=len(null_win))
    factors = []
    for name, pl in pools.items():
        fit_pool(pl, tx, region, zone, bins, bg, null_pos, null_win, n_tx)
        factors.append(pd.DataFrame({"pool": name, "bin": bin_names(), "factor": pl.f}))
        stats.update({f"prior_shape {name} {r}": a for r, a in zip(REGIONS, pl.a)})
    null = pooled_moments(pools, null_win, tx, region, zone, bins)
    null_region = region[null_win]
    phi = calibrate_phi(null.h, null.mean, null.var, null_region)
    stats.update(phi_stats(phi))

    # the test of the candidates; BH per class, over every candidate of the class
    if scan:
        pos = scan_candidates(layout, tx, pools, cs_pos, start_codons)
        stats.update(start_codons=",".join(start_codons), candidates=len(pos))
    else:
        pos, ids, got = given_starts(orfs_path, layout)
        stats.update(got)
    h, e, m, var = pooled_moments(pools, pos, tx, region, zone, bins)
    p = pvalues(h, m, var, region[pos], phi)
    cls, n_tests, tested = class_tests(layout, tx, region, pos, cs_pos, scan, start_codons)
    if scan:
        # null windows clear of the called starts: harringtonine raises their codons 1-10, as for the annotated
        # starts. Call, drop the null windows in RAISED of the called starts, calibrate phi again, call again
        first_calls, _ = near_starts(pos, p, class_q(p, cls, n_tests) < fdr, np.flatnonzero(cls == ANNOTATED))
        keep = ~covered(layout.size, pos[first_calls] + RAISED[0], pos[first_calls] + RAISED[1])[null_win]
        null, null_region = Moments(*(x[keep] for x in null)), null_region[keep]
        phi = calibrate_phi(null.h, null.mean, null.var, null_region)
        p = pvalues(h, m, var, region[pos], phi)
        stats.update(phi_stats(phi))
        stats["null_windows_near_called"] = int((~keep).sum())
    null_p = pvalues(null.h, null.mean, null.var, null_region, phi)
    ann = np.isin(pos, cs_pos) & (m >= 3)
    enrich = np.median(h[ann] / m[ann]) if ann.any() else np.nan
    stats["typical_enrichment"] = enrich
    typical_p = (pvalues(np.ceil(enrich * m).astype(np.int64), m, var, region[pos], phi) if ann.any()
                 else np.full(len(pos), np.nan))
    t = tx[pos].astype(np.int64)
    start = pos - layout.first[t]
    out = pd.DataFrame({
        "Name": layout.names[t], "start": start,
        "codon": [layout.seq[x:x + 3].tobytes().decode() for x in pos],
        "region": np.array(REGIONS)[region[pos]],
        "frame": pd.Series((start - layout.cds_start[t]) % 3, dtype="Int64").mask(layout.cds_start[t] < 0),
        "ORF_id": [f"{n}_{x + 1}" for n, x in zip(layout.names[t], start)] if scan else ids,
        "harringtonine": h, "elongation": e, "expected": m,
        "enrichment": h / np.maximum(m, 1e-300), "p": p, "typical_p": typical_p,
        "multimapped": sum(window_sums(pl.multi, pos) for pl in pools.values())
        / np.maximum(sum(window_sums(pl.h, pos) + window_sums(pl.e, pos) for pl in pools.values()), 1),
        "shoulder": shoulder(pools.values(), pos),
    })
    q = class_q(p, cls, n_tests)
    out["q"] = q
    stats.update(class_stats(p, q, cls, n_tests, null_p, null_region, fdr))
    columns = STARTS_COLUMNS
    if scan:
        orfs, catalogue = scan_calls(out, pos, p, q, cls, layout, tx, fdr, stats)
        write_tsv(catalogue, f"{out_prefix}.catalogue.tsv", na_rep="NA")
        write_tsv(orfs, f"{out_prefix}.start_orfs.tsv", na_rep="NA")
        columns = STARTS_COLUMNS + SCAN_COLUMNS
    called = out["called"].to_numpy() if scan else q < fdr
    write_tsv(codon_table(layout, tested, pos, cls, called), f"{out_prefix}.start_codons.tsv")
    out = out[columns]
    log.info("%s", stats)
    write_tsv(out, f"{out_prefix}.starts.tsv", na_rep="NA")
    write_tsv(kernels, f"{out_prefix}.kernel.tsv")
    write_tsv(pd.concat(factors, ignore_index=True), f"{out_prefix}.factors.tsv")
    finish(stats, out_prefix, report)
    return out, stats

import filecmp

import numpy as np
import pandas as pd
import pysam
import pytest
from scipy import stats as sps

from pitstop import cli, orfs, starts
from pitstop.quant import quantify

from synth import make_start_dataset

WINDOW = (28, 30)
PLANTED = ["S1_uATG", "S1_uCTG", "S2_uoORF", "S3_ext", "S4_pairATG", "S5_first", "N1_ORF"]


def test_nb_sf_matches_scipy():
    rng = np.random.default_rng(0)
    h = rng.integers(0, 30, 500)
    mean = rng.random(500) * 20
    size = rng.random(500) * 5 + 0.1
    size[:100] = np.inf
    fin = np.isfinite(size)
    ref = np.where(fin, sps.nbinom.sf(h - 1, np.where(fin, size, 1), np.where(fin, size, 1) / (np.where(fin, size, 1) + mean)),
                   sps.poisson.sf(h - 1, mean))
    assert np.allclose(starts.nb_sf(h, mean, size), ref, atol=1e-12)


def test_null_codons_keep_clear_of_starts():
    bg = np.zeros(200, dtype=bool)
    bg[64:150] = True
    is_start = np.zeros(200, dtype=bool)
    is_start[100] = True
    got = set(starts.null_codons(bg, is_start).tolist())
    # window [c - 1, c + 3) inside the background, no start within 8 nt
    assert got == {c for c in range(65, 148) if abs(c - 100) > 8}


def test_fit_factors_recover_planted():
    rng = np.random.default_rng(1)
    n_tx, per = 400, 60
    tx = np.repeat(np.arange(n_tx), per)
    b = np.tile(np.repeat(np.arange(3), per // 3), n_tx)
    lam = rng.gamma(2, 1, tx.size)
    f_true, s_true = np.array([4.0, 1.0, 0.5]), rng.uniform(0.5, 2, n_tx)
    e = rng.poisson(lam)
    h = rng.poisson(s_true[tx] * f_true[b] * lam)
    pos = np.arange(tx.size)
    f, s = starts.fit_factors(h, e, tx, b, pos, pos, n_tx, 3, prior_reads=0)
    assert f / f[1] == pytest.approx(f_true, rel=0.1)


def test_calibrate():
    # per read length: P-sites at nt -15..62 of the annotated starts (a peak on a flat floor) and CDS-interior P-sites
    # per frame. The start peak sets the shift; the frame only vetoes it when it is clearly elsewhere
    def prof(peak, at=100):
        r = np.ones(starts.PROFILE[1] - starts.PROFILE[0])
        r[peak - starts.PROFILE[0]] = at
        return r
    cases = [  # (start peak, P-sites at it, frames) -> (calibrated, shift)
        (0, 100, [900, 50, 50], True, 0),
        (-1, 100, [300, 100, 600], True, 1),       # 5' ends 1 nt longer: frame 2 agrees
        (2, 100, [100, 100, 800], True, -2),       # frame 2 too, but a codon later: the frame alone would say +1
        (-3, 100, [900, 50, 50], False, 0),        # too far; the clear frame 0 keeps it in place
        (0, 10, [100, 50, 850], False, 1),         # too few P-sites at the peak; the clear frame 2 moves it
        (-1, 100, [800, 100, 100], False, 0),      # the frame is clearly against the peak
        (-1, 100, [480, 40, 480], True, 1),        # frames 0 and 2 alike: no veto
        (-1, 100, [5, 1, 60], True, 1)]            # too few interior P-sites for the frame to decide
    c = starts.calibrate(np.array([prof(p, a) for p, a, *_ in cases]), np.array([f for _, _, f, *_ in cases]),
                         (24, 24 + len(cases) - 1))
    assert c["length"].tolist() == list(range(24, 24 + len(cases)))
    assert c["calibrated"].tolist() == [x[3] for x in cases]
    assert c["shift"].tolist() == [x[4] for x in cases]
    assert not c["mixed"].any()


def test_calibrate_mixed():
    # a nt beside the start peak with half the peak's P-sites or more: the 5' ends split between two offsets, the
    # length is mixed (it still moves by minus the peak)
    def prof(peak, at, beside, side):
        r = np.ones(starts.PROFILE[1] - starts.PROFILE[0])
        r[peak - starts.PROFILE[0]] = at
        r[peak + side - starts.PROFILE[0]] = beside
        return r
    cases = [  # (start peak, P-sites at it, P-sites beside it, side) -> mixed
        (-1, 100, 88, 1, True),        # eIF 30 nt: offsets 13 and 12
        (0, 100, 45, -1, False),       # eIF 25 nt at most
        (0, 100, 50, -1, True),
        (-15, 100, 60, 1, True),       # a peak at the profile's first nt
        (0, 10, 9, 1, False)]          # too few P-sites at the peak to tell
    c = starts.calibrate(np.array([prof(*x[:4]) for x in cases]), np.ones((len(cases), 3)),
                         (24, 24 + len(cases) - 1))
    assert c["start_peak"].tolist() == [x[0] for x in cases]
    assert c["beside_peak"].tolist() == [x[2] for x in cases]
    assert c["beside_nt"].tolist() == [x[0] + x[3] for x in cases]
    assert np.allclose(c["beside_share"], [x[2] / x[1] for x in cases])
    assert c["mixed"].tolist() == [x[4] for x in cases]
    assert not c["pair"].any()
    assert c.loc[0, "calibrated"] and c.loc[0, "shift"] == 1


def test_calibrate_pair():
    # a pair length: the anchor is the start peak and the larger nt beside it, and the pair moves onto nt 0 and 1,
    # whichever of the two is larger; the frame vetoes only when it is in neither frame of the pair
    def prof(peak, at, beside, side):
        r = np.ones(starts.PROFILE[1] - starts.PROFILE[0])
        r[peak - starts.PROFILE[0]] = at
        r[peak + side - starts.PROFILE[0]] = beside
        return r
    cases = [  # (start peak, P-sites at it, P-sites beside it, side, frames) -> (calibrated, shift)
        (-1, 100, 88, 1, [530, 80, 390], True, 1),     # eIF 30 nt: the pair at nt -1 and 0
        (0, 100, 88, -1, [530, 80, 390], True, 1),     # the same pair, the other nt larger: the same shift
        (0, 100, 70, 1, [360, 550, 90], True, 0),      # K562 27 nt: frame 1 leads, but it is a frame of the pair
        (1, 100, 80, -1, [240, 720, 40], True, 0),     # the same, the peak at +1
        (2, 100, 60, -1, [300, 100, 600], True, -1),   # the pair at 1 and 2
        (-2, 100, 60, 1, [300, 100, 600], True, 2),    # the pair at -2 and -1: the farthest it moves
        (-3, 100, 60, 1, [500, 50, 450], False, 1),    # too far; the clear frame pair (2, 0) moves it
        (3, 100, 60, -1, [900, 50, 50], False, 0),     # a codon-1 peak: too far, though its pair starts at +2
        (-1, 100, 80, 1, [50, 900, 50], False, 0),     # frame 1 is in neither frame of the pair; no clear frame pair
        (0, 10, 9, 1, [5, 1, 60], False, 0)]           # too few P-sites at the peak and inside the CDSs
    p = np.array([prof(*x[:4]) for x in cases])
    f = np.array([x[4] for x in cases])
    lengths = (24, 24 + len(cases) - 1)
    c = starts.calibrate(p, f, lengths, pair=np.ones(len(cases), dtype=bool))
    assert c["pair"].all()
    assert c["calibrated"].tolist() == [x[5] for x in cases]
    assert c["shift"].tolist() == [x[6] for x in cases]
    single = starts.calibrate(p, f, lengths)   # the peak alone: the coin flip of the first two, K562's veto
    assert single["shift"].tolist()[:3] == [1, 0, -1] and not single.loc[2, "calibrated"]


def test_pair_lengths():
    # a pair length is mixed in half of the harringtonine libraries or more
    rows = pd.DataFrame({"length": [27, 27, 27, 28, 28, 28, 29, 29, 29],
                         "mixed": [True, True, False, True, False, False, False, False, False]})
    assert starts.pair_lengths(rows).tolist() == [True, False, False]
    assert starts.pair_lengths(rows[rows.index % 3 != 2]).tolist() == [True, True, False]   # 1 of 2


def test_choose_start_lengths():
    # two libraries per length (a harringtonine and an elongation one): the window, calibration in the harringtonine
    # library and a frame in both decide; of the passing runs, the one with the most P-sites. Mixed offsets do not
    rows = []
    for length, cal, mix, per, reads in [(22, True, "", True, 900), (23, True, "harringtonine", True, 10),
                                         (24, False, "", True, 10), (25, True, "elongation", True, 25),
                                         (26, True, "", True, 25), (27, True, "", False, 10), (28, True, "", True, 20),
                                         (29, True, "", True, 20), (30, True, "", True, 20), (35, True, "", True, 10)]:
        for kind in ("harringtonine", "elongation"):
            rows.append({"library": kind, "kind": kind, "length": length, "psites": reads,
                         "calibrated": cal or kind == "elongation", "mixed": kind == mix, "periodic": per})
    why, best = starts.choose_start_lengths(pd.DataFrame(rows), (23, 34))
    assert best == [28, 29, 30]          # 60 P-sites per library against 25-26's 50
    assert why[22] == why[35] == "outside the window 23-34"
    assert why[23] == "outside the run 28-30"
    assert why[24] == "not calibrated in 1 of 1 harringtonine libraries"
    assert why[27] == "no clear frame in 2 of 2 libraries"
    assert why[25] == why[26] == "outside the run 28-30" and why[28] == "pass"   # elongation libraries: no say
    assert starts.choose_start_lengths(pd.DataFrame(rows), (31, 34))[1] == []


def test_shift_psites_stay_in_their_transcript():
    tx = np.array([-1, 0, 0, 0, 0, -1, 1, 1, 1, -1])         # two transcripts and the gaps around them
    pos, length = np.array([1, 2, 4, 7, 8]), np.array([28, 28, 29, 28, 29])
    multi = np.array([False, True, False, True, False])
    got = starts.shift_psites((pos, length, multi), np.array([-1, 1]), 28, tx)   # 28 nt: -1, 29 nt: +1
    assert [g.tolist() for g in got] == [[1, 6], [28, 28], [True, True]]


def test_orf_types():
    cs, ce = np.full(8, 100), np.full(8, 400)
    start = np.array([10, 40, 97, 98, 100, 160, 161, 420])
    end = np.array([60, 400, 400, 130, 400, 400, 200, 500])
    assert list(starts.orf_types(start, end, cs, ce)) == [
        "uORF", "extension", "extension", "uoORF", "CDS", "truncation", "internal", "dORF"]
    assert starts.orf_types(np.array([5]), np.array([50]), np.array([-1]), np.array([-1]))[0] == "other"


def test_bh_counts_untested_as_p1():
    p = np.array([0.001, 0.02, 0.5])
    assert np.allclose(starts.bh(p, 3), sps.false_discovery_control(p))
    assert np.allclose(starts.bh(p, 10), sps.false_discovery_control(np.r_[p, np.ones(7)])[:3])


def test_near_starts_best_first():
    pos = np.array([100, 106, 112, 150, 200])
    p = np.array([1e-3, 1e-8, 1e-5, 1e-6, 0.5])
    called, by = starts.near_starts(pos, p, p < 0.01)
    # 106 is best: it stops 100 (6 nt before) and 112 (codon 2 after); 150 is 44 nt after, out of reach
    assert called.tolist() == [False, True, False, True, False]
    assert by.tolist() == [1, -1, 1, -1, -1]
    # codons 1-10 of an annotated start (here 100) are not called, even when it is not called itself
    called, by = starts.near_starts(pos, p, p < 0.01, annotated=[0])
    assert called.tolist() == [True, False, False, True, False] and by.tolist() == [-1, 0, 0, -1, -1]


def test_shoulder():
    # log2 (h + 1) / (e + 1) in the codons 1-5 after a start (nt 3-17), minus the same in nt -18..-4; the reads of the
    # pools add before the ratio. The window and the nt just outside the two zones do not count
    s = 50
    h, e = np.zeros(100, dtype=np.int64), np.zeros(100, dtype=np.int64)
    h[[s + 3, s + 17]], e[s + 10] = 3, 1
    h[s - 18], e[s - 4] = 1, 3
    out = [s - 19, s - 3, s - 1, s, s + 2, s + 18]
    h[out], e[out] = 100, 50
    pools = [starts.Pool(h, e, None), starts.Pool(h, e, None)]
    assert np.allclose(starts.shoulder(pools, np.array([s])), np.log2(13 / 3) - np.log2(3 / 7))


def pool_sheet(d, data):
    """`orfs` on each library with one offset table per pool, from `quant` on the merged elongation BAMs (as
    docs/starts.md advises); the library sheet of `starts`."""
    pysam.cat("-o", str(d / "elong.bam"), *[str(data["libs"]["elongation", r]) for r in (1, 2)])
    quantify(str(d / "elong.bam"), str(data["gtf"]), str(data["fasta"]), WINDOW, str(d / "pool"))
    rows = []
    for rep in (1, 2):
        for kind in ("harringtonine", "elongation"):
            orfs.run(str(data["libs"][kind, rep]), str(data["gtf"]), str(data["fasta"]), WINDOW, str(d / f"{kind}{rep}"),
                     orfs_path=str(data["orfs"]), offsets_path=str(d / "pool.offsets.tsv"))
        rows.append(("A", str(d / f"harringtonine{rep}"), str(d / f"elongation{rep}")))
    sheet = d / "sheet.tsv"
    pd.DataFrame(rows, columns=["pool", "harringtonine", "elongation"]).to_csv(sheet, sep="\t", index=False)
    return sheet


def orfs_offsets(qc):
    """Per row of lengths.tsv: the offset of its length in its library's orfs run (frames.tsv)."""
    return np.array([pd.read_csv(f"{lib}.frames.tsv", sep="\t").set_index("length").at[length, "offset"]
                     for lib, length in zip(qc["library"], qc["length"])])


@pytest.fixture(scope="module")
def started(tmp_path_factory):
    d = tmp_path_factory.mktemp("startsynth")
    data = make_start_dataset(d)
    sheet = pool_sheet(d, data)
    out, st = starts.run(str(sheet), str(data["gtf"]), str(data["fasta"]), str(d / "s"), str(data["orfs"]))
    scan, scan_st = starts.run(str(sheet), str(data["gtf"]), str(data["fasta"]), str(d / "scan"), scan=True)
    return dict(dir=d, data=data, sheet=sheet, out=out, stats=st, scan=scan, scan_stats=scan_st)


def test_start_lengths_from_kernel(started):
    assert started["stats"]["start_lengths"] == "28,29,30"
    k = pd.read_csv(started["dir"] / "s.kernel.tsv", sep="\t")
    peak = k.loc[k.groupby("length")["kernel"].idxmax(), "nt"]
    assert (peak == 0).all()
    qc = pd.read_csv(started["dir"] / "s.lengths.tsv", sep="\t")
    assert list(qc.columns) == starts.LENGTH_COLUMNS and len(qc) == 4 * 3
    assert qc["calibrated"].all() and (qc["shift"] == 0).all() and qc["periodic"].all() and qc["start_length"].all()
    assert not qc["mixed"].any() and (qc["beside_peak"] < starts.MIXED * qc["at_peak"]).all()
    assert (qc["offset"] == orfs_offsets(qc)).all()          # no shift: the orfs runs' offsets
    assert (qc["qc"] == "pass").all() and (qc["pool"] == "A").all()
    assert (qc["psites"] == qc[["5UTR", "CDS", "3UTR", "transcript"]].sum(axis=1)).all()
    assert np.allclose(qc.groupby("library")["share"].sum(), 1)
    assert np.allclose(qc.groupby("library")["interior_share"].sum(), 1)
    prof = pd.read_csv(started["dir"] / "s.profiles.tsv", sep="\t")
    n = (starts.PROFILE[1] - starts.PROFILE[0]) + (starts.STOP_PROFILE[1] - starts.STOP_PROFILE[0]) + starts.CDS_BINS * 3
    assert len(prof) == 4 * 3 * n
    start = prof[prof["anchor"] == "start"].merge(qc, on=["library", "length"])
    assert (start.loc[start["x"] == start["start_peak"], "psites_x"] == start.loc[start["x"] == start["start_peak"],
                                                                                   "at_peak"]).all()
    cds = prof[prof["anchor"] == "cds"].groupby(["library", "length"])["psites"].sum()
    assert (cds - qc.set_index(["library", "length"])["CDS"]).eq(0).all()   # aligned on (library, length)


@pytest.fixture(scope="module")
def extended(tmp_path_factory):
    # the harringtonine libraries leave most 30 nt 5' ends 1 nt longer (another digestion): with the pool's offsets
    # (from the elongation libraries) their P-sites peak at nt -1 and sit in frame 2 inside the CDSs
    d = tmp_path_factory.mktemp("extended")
    data = make_start_dataset(d, harr_offsets={30: [(13, .7), (12, .2), (14, .1)]})
    sheet = pool_sheet(d, data)
    out, st = starts.run(str(sheet), str(data["gtf"]), str(data["fasta"]), str(d / "s"), str(data["orfs"]))
    return dict(dir=d, out=out, stats=st)


def test_calibration_moves_a_shifted_length(extended):
    # the 30 nt P-sites of the harringtonine libraries move by +1, nothing else moves; 30 nt stays a start length
    # (the peak-at-nt-0 rule kept 28-29 only) and the calls are those of an unshifted set
    cal = pd.read_csv(extended["dir"] / "s.lengths.tsv", sep="\t")
    moved = (cal["kind"] == "harringtonine") & (cal["length"] == 30)
    assert (cal.loc[moved, "start_peak"] == -1).all() and (cal.loc[moved, "shift"] == 1).all()
    assert (cal.loc[moved, "frame2"] > cal.loc[moved, ["frame0", "frame1"]].max(axis=1)).all()
    assert (cal.loc[~moved, "shift"] == 0).all() and cal["calibrated"].all()
    assert (cal["offset"] == orfs_offsets(cal) + moved).all()   # 30 nt of the harringtonine libraries: 1 nt more
    assert not cal["mixed"].any()                  # 20% of the 5' ends at offset 12: one offset holds
    assert extended["stats"]["start_lengths"] == "28,29,30"
    k = pd.read_csv(extended["dir"] / "s.kernel.tsv", sep="\t")
    assert (k.loc[k.groupby("length")["kernel"].idxmax(), "nt"] == 0).all()
    out = extended["out"].set_index("ORF_id")
    q = out["q"]
    assert (out.loc[out.index.str.endswith(":CDS"), "q"] < 0.05).all()
    assert sum(q[i] < 0.05 for i in PLANTED) >= 6
    assert q["S6_pause"] > 0.05 and q["S7_quiet"] > 0.05


def test_mixed_length_anchors_its_pair(tmp_path):
    # the harringtonine libraries split their 30 nt 5' ends between offsets 13 and 12 (as eIF): the start peak has
    # the nt beside it at about its height, so 30 nt is mixed there (the elongation libraries' is not: 0.2 of its 5'
    # ends at 13) and a pair length in every library. The pair (nt -1 and 0, whichever is larger) moves onto nt 0 and
    # 1 in all four libraries, and 30 nt is a start length; the calls are those of an unshifted set
    data = make_start_dataset(tmp_path, harr_offsets={30: [(13, .5), (12, .5)]})
    sheet = pool_sheet(tmp_path, data)
    out, st = starts.run(str(sheet), str(data["gtf"]), str(data["fasta"]), str(tmp_path / "s"), str(data["orfs"]))
    qc = pd.read_csv(tmp_path / "s.lengths.tsv", sep="\t")
    thirty = qc["length"] == 30
    mixed = (qc["kind"] == "harringtonine") & thirty
    assert qc.loc[mixed, "mixed"].all() and not qc.loc[~mixed, "mixed"].any()
    assert qc.loc[thirty, "pair"].all() and not qc.loc[~thirty, "pair"].any()
    assert qc.loc[mixed, "start_peak"].isin([-1, 0]).all() and qc.loc[mixed, "beside_nt"].isin([-1, 0]).all()
    assert qc["calibrated"].all() and (qc.loc[thirty, "shift"] == 1).all() and (qc.loc[~thirty, "shift"] == 0).all()
    assert (qc["offset"] == orfs_offsets(qc) + thirty).all()
    assert (qc["qc"] == "pass").all() and st["start_lengths"] == "28,29,30"
    out = out.set_index("ORF_id")
    q = out["q"]
    assert (out.loc[out.index.str.endswith(":CDS"), "q"] < 0.05).all()
    assert sum(q[i] < 0.05 for i in PLANTED) >= 6
    assert q["S6_pause"] > 0.05 and q["S7_quiet"] > 0.05
    # the report: the pair's move is not an offset of the initiation library
    doc = (tmp_path / "s.starts.html").read_text(encoding="utf-8")
    assert "Mixed lengths: 30 nt. All are pair lengths." in doc and "30 nt moved +1 (pair)" in doc
    assert "The calibration moves the pair lengths 30 nt by +1 in the initiation library" in doc
    assert "Some initiation libraries need other offsets." not in doc


def test_synthetic_calls(started):
    # synthetic truth for given starts. Over seeds 3-12: every annotated start called; 1 of 70 planted starts missed
    # (a uORF CTG, q 0.99); the untranslated uORF never called, the pause in 3 seeds (q 0.015-0.040); at most two
    # of the ten untranslated filler candidates called (q 0.027-0.043: BH among true starts). The pair's CTG gets
    # the ATG's nt -1 reads but is never called (the scan's near-start rule would drop it). The pair's ATG and the
    # first of the two close starts always have the smaller p.
    out = started["out"].set_index("ORF_id")
    q, p = out["q"], out["p"]
    cds = out[out.index.str.endswith(":CDS")]
    assert (cds["q"] < 0.05).all()
    assert sum(q[i] < 0.05 for i in PLANTED) >= 6
    assert q["S6_pause"] > 0.05 and q["S7_quiet"] > 0.05
    assert (q[out.index.str.endswith("_cand")] < 0.05).sum() <= 1
    assert p["S4_pairATG"] < p["S4_pairCTG"] and p["S5_first"] < p["S5_after"]
    assert out.at["S6_pause", "enrichment"] < 1.5        # a pause: elongation peak, harringtonine not raised
    assert out.at["S2_uoORF", "region"] == "5UTR" and pd.isna(out.at["N1_ORF", "frame"])


def test_starts_cli_deterministic(started, tmp_path):
    d = started["dir"]
    data = started["data"]
    cli.main(["starts", "--libraries", str(started["sheet"]), "--gtf", str(data["gtf"]), "--fasta", str(data["fasta"]),
              "--orfs", str(data["orfs"]), "--out-prefix", str(tmp_path / "again")])
    for f in ("starts", "lengths", "profiles", "libraries", "pairs", "start_codons", "kernel", "factors"):
        assert filecmp.cmp(d / f"s.{f}.tsv", tmp_path / f"again.{f}.tsv", shallow=False), f
    stats = [pd.read_csv(x, sep="\t", dtype=str).set_index("stat")["value"] for x in (d / "s.starts_stats.tsv",
                                                                                     tmp_path / "again.starts_stats.tsv")]
    assert stats[0].equals(stats[1])
    cli.main(["starts", "--libraries", str(started["sheet"]), "--gtf", str(data["gtf"]), "--fasta", str(data["fasta"]),
              "--scan", "--out-prefix", str(tmp_path / "scan")])
    for f in ("starts", "start_orfs", "start_codons", "starts_stats"):
        assert filecmp.cmp(d / f"scan.{f}.tsv", tmp_path / f"scan.{f}.tsv", shallow=False), f


def test_pairs_and_libraries(started):
    # a pair per row of the sheet, its distance from the two libraries' interior shares; harringtonine raises the start
    # peak (synthetic: the start profiles of harringtonine have 5x the elongation's peak)
    d = started["dir"]
    pairs = pd.read_csv(d / "s.pairs.tsv", sep="\t")
    sheet = pd.read_csv(started["sheet"], sep="\t")
    assert list(pairs.columns) == starts.PAIR_COLUMNS and (pairs[["pool", "harringtonine", "elongation"]]
                                                           == sheet).all().all()
    share = pd.read_csv(d / "s.lengths.tsv", sep="\t").pivot(index="length", columns="library", values="interior_share")
    for _, r in pairs.iterrows():
        assert r["distance"] == pytest.approx(0.5 * (share[r["harringtonine"]] - share[r["elongation"]]).abs().sum())
    assert (pairs["distance"] < starts.PAIR_DISTANCE).all()          # one digestion
    libs = pd.read_csv(d / "s.libraries.tsv", sep="\t").set_index("library")
    assert list(libs.columns) == starts.LIBRARY_COLUMNS[1:]
    assert np.allclose(libs[["5UTR", "CDS", "3UTR", "transcript"]].sum(axis=1), 1)
    assert np.allclose(libs[["frame0", "frame1", "frame2"]].sum(axis=1), 1)
    for _, r in pairs.iterrows():
        assert libs.at[r["harringtonine"], "start_enrichment"] > 3 * libs.at[r["elongation"], "start_enrichment"]


def test_start_codons(started):
    # per class: the tests are the class's tests in starts_stats.tsv, the calls the scan's calls; shares sum to 1
    for name, st in (("s", started["stats"]), ("scan", started["scan_stats"])):
        c = pd.read_csv(started["dir"] / f"{name}.start_codons.tsv", sep="\t", keep_default_na=False)
        assert list(c.columns) == starts.CODON_COLUMNS
        g = c.groupby("class")
        for cls, tests in g["tests"].sum().items():
            assert tests == st[f"tests {cls}"], (name, cls)
        for what in ("test_share", "call_share"):
            sums = g[what].sum()
            assert np.allclose(sums[sums > 0], 1)
    c = pd.read_csv(started["dir"] / "scan.start_codons.tsv", sep="\t", keep_default_na=False)
    assert c["called"].sum() == started["scan_stats"]["called"]
    assert set(c["codon"]) <= set(starts.START_CODONS)
    ann = c[c["class"] == "annotated"].set_index("codon")
    assert ann.at["ATG", "test_share"] == 1             # synthetic annotated starts are ATG


def test_length_window(started, tmp_path):
    # a window that holds no read length of the data: no start lengths, but lengths.tsv says why; --start-lengths
    # still overrides the QC
    data = started["data"]
    args = ["starts", "--libraries", str(started["sheet"]), "--gtf", str(data["gtf"]), "--fasta", str(data["fasta"]),
            "--orfs", str(data["orfs"])]
    with pytest.raises(ValueError, match="no read length passes the QC"):
        cli.main(args + ["--out-prefix", str(tmp_path / "w"), "--length-window", "31-34"])
    qc = pd.read_csv(tmp_path / "w.lengths.tsv", sep="\t")
    assert qc["qc"].str.startswith("outside the window 31-34").all() and not qc["start_length"].any()
    cli.main(args + ["--out-prefix", str(tmp_path / "g"), "--length-window", "31-34", "--start-lengths", "29-30"])
    qc = pd.read_csv(tmp_path / "g.lengths.tsv", sep="\t")
    assert qc.loc[qc["start_length"], "length"].unique().tolist() == [29, 30]
    assert pd.read_csv(tmp_path / "g.starts_stats.tsv", sep="\t").set_index("stat").at["start_lengths", "value"] == "29,30"


def test_scan_calls(started):
    # synthetic truth for the scan. Over seeds 3-12: every annotated start called; 66 of 70 planted starts called at their
    # codon, the misses at 1.1-4.2x enrichment with 3-6 elongation reads (q 0.08-1); the pair's CTG and the weaker
    # start 5 codons after a planted one are never called, and listed with the start that stops them; the pause
    # and the untranslated uORF never called. Calls at no planted start: 0-4 per seed in 5' UTRs (q >= 0.012 but
    # one), 1-6 in CDSs
    data = started["data"]
    out = started["scan"]
    key = {f"{n}_{s + 1}": i for i, (n, s, e, _, _) in data["orf_coords"].items()}
    out = out.assign(truth=out["ORF_id"].map(key)).set_index("ORF_id")
    by_truth = out.dropna(subset=["truth"]).set_index("truth")
    cds = by_truth[by_truth.index.str.endswith(":CDS")]
    assert cds["called"].all() and (cds["type"] == "CDS").all()
    assert sum(by_truth.at[i, "called"] for i in PLANTED) >= 6
    want = {"S1_uATG": "uORF", "S1_uCTG": "uORF", "S2_uoORF": "uoORF", "S3_ext": "extension", "S4_pairATG": "uORF",
            "S5_first": "uORF", "N1_ORF": "other"}
    for i, t in want.items():
        assert by_truth.at[i, "type"] == t, i
        assert by_truth.at[i, "end"] == data["orf_coords"][i][2], i
    for weaker, stronger in (("S4_pairCTG", "S4_pairATG"), ("S5_after", "S5_first")):
        if weaker in by_truth.index:
            assert not by_truth.at[weaker, "called"]
            if by_truth.at[stronger, "called"]:
                assert out.at[by_truth.at[weaker, "stopped_by"], "truth"] == stronger
    for i in ("S6_pause", "S7_quiet"):
        assert i not in by_truth.index or not by_truth.at[i, "called"], i
    false_5utr = out[out["called"] & out["truth"].isna() & (out["region"] == "5UTR")]
    assert len(false_5utr) <= 3
    so = pd.read_csv(started["dir"] / "scan.start_orfs.tsv", sep="\t")
    assert len(so) == started["scan_stats"]["called"] - len(cds) - started["scan_stats"]["called_without_stop"]
    assert (so["type"] != "CDS").all() and set(so.columns[:4]) == {"ORF_id", "Name", "start", "end"}


def test_scan_start_codons_and_fdr(started, tmp_path):
    # --start-codons chooses the scan's candidates and BH counts only them; the annotated starts stay. --fdr is the
    # cut of a call
    data = started["data"]
    cli.main(["starts", "--libraries", str(started["sheet"]), "--gtf", str(data["gtf"]), "--fasta", str(data["fasta"]),
              "--scan", "--start-codons", "atg", "--fdr", "0.01", "--out-prefix", str(tmp_path / "atg")])
    out = pd.read_csv(tmp_path / "atg.starts.tsv", sep="\t")
    st = pd.read_csv(tmp_path / "atg.starts_stats.tsv", sep="\t").set_index("stat")["value"]
    cds = out["type"] == "CDS"
    assert (out.loc[~cds, "codon"] == "ATG").all() and cds.sum() == (started["scan"]["type"] == "CDS").sum()
    assert (out.loc[out["called"], "q"] < 0.01).all() and st["start_codons"] == "ATG" and "q<0.01 5UTR" in st.index
    assert 0 < int(st["tests 5UTR"]) < started["scan_stats"]["tests 5UTR"]
    key = {f"{n}_{s + 1}": i for i, (n, s, e, _, _) in data["orf_coords"].items()}
    called = set(out.loc[out["called"], "ORF_id"].map(key).dropna())
    assert {"S1_uATG", "S4_pairATG"} <= called and "S1_uCTG" not in called
    with pytest.raises(ValueError, match="scan"):
        starts.run(str(started["sheet"]), str(data["gtf"]), str(data["fasta"]), str(tmp_path / "x"),
                   str(data["orfs"]), start_codons=("ATG",))


def test_scan_null_windows_clear_of_called_starts(started, tmp_path, monkeypatch):
    # a scan calibrates phi again on the null windows outside the codons 1-10 of the starts of a first call (run-off,
    # as after an annotated start); given starts calibrate it once. The synthetic sets have too few null windows near
    # a start for phi to change (seeds 3-12: 0-6 of 281-353), so a loose cut makes the first call reach them
    seen = []
    real = starts.calibrate_phi
    monkeypatch.setattr(starts, "calibrate_phi", lambda h, m, v, r: seen.append(len(h)) or real(h, m, v, r))
    data = started["data"]
    _, st = starts.run(str(started["sheet"]), str(data["gtf"]), str(data["fasta"]), str(tmp_path / "loose"), scan=True,
                       fdr=0.5)
    assert st["null_windows_near_called"] > 0
    assert seen == [st["null_windows"], st["null_windows"] - st["null_windows_near_called"]]
    seen.clear()
    _, st = starts.run(str(started["sheet"]), str(data["gtf"]), str(data["fasta"]), str(tmp_path / "given"),
                       str(data["orfs"]))
    assert seen == [st["null_windows"]] and "null_windows_near_called" not in st


@pytest.mark.parametrize("table", ["start_orfs", "catalogue"])
def test_scan_orfs_feed_orfs(started, tmp_path, table):
    # start_orfs.tsv and catalogue.tsv are ORF tables for `pitstop orfs`: every row is kept, and the catalogue's CDS
    # rows are the annotated CDSs (none added)
    d, data = started["dir"], started["data"]
    _, _, _, st = orfs.run(str(data["libs"]["elongation", 1]), str(data["gtf"]), str(data["fasta"]), WINDOW,
                           str(tmp_path / "o"), orfs_path=str(d / f"scan.{table}.tsv"),
                           offsets_path=str(d / "pool.offsets.tsv"))
    assert st["orf_table_rows"] == len(pd.read_csv(d / f"scan.{table}.tsv", sep="\t"))
    assert sum(v for k, v in st.items() if k.startswith("orfs_dropped_")) == 0
    if table == "catalogue":
        assert st["orfs_cds_added"] == 0


def test_catalogue(started):
    out = started["scan"]
    cat = pd.read_csv(started["dir"] / "scan.catalogue.tsv", sep="\t").set_index("ORF_id")
    calls = out[out["called"] & out["type"].isin(["uORF", "uoORF"])]
    assert set(cat.index[cat["type"] != "CDS"]) == set(calls["ORF_id"])
    assert (cat["type"] == "CDS").sum() == 37
    # the shoulder of each row is its start's in starts.tsv (synthetic 5' UTRs have no shoulder contrast: their
    # harringtonine / elongation ratio, HARR_5UTR, is about the planted kernel's in codons 1-5)
    shoulder = out.set_index("ORF_id")["shoulder"]
    assert np.allclose(cat["shoulder"], shoulder[cat.index]) and shoulder.notna().all()
    assert started["out"]["shoulder"].notna().all()
    data = started["data"]
    for i in ("S1_uATG", "S2_uoORF", "N1_ORF", "S3_ext"):
        n, s, e, _, _ = data["orf_coords"][i]
        assert (f"{n}_{s + 1}" in cat.index) == (i in ("S1_uATG", "S2_uoORF")), i   # no extension, no ORF without a CDS

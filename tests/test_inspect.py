import filecmp

import numpy as np
import pandas as pd
import pytest

from pitstop import annotation, bam, inspect

from synth import make_inspect_dataset

TABLES = ("lengths", "anchors", "ends", "complexity", "summary")


def run_inspect(data, prefix):
    s = inspect.run(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), str(prefix), report=False)
    t = {name: pd.read_csv(f"{prefix}.{name}.tsv", sep="\t", dtype={"length": str} if name == "ends" else None)
         for name in TABLES}
    t["lengths"] = t["lengths"].set_index("length")
    t["table"] = t.pop("summary")
    return dict(prefix=prefix, summary=s, **t)


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
    data = make_inspect_dataset(tmp_path_factory.mktemp("clean"))
    return data, run_inspect(data, tmp_path_factory.mktemp("clean_run") / "lib")


@pytest.fixture(scope="module")
def untemplated(tmp_path_factory):
    data = make_inspect_dataset(tmp_path_factory.mktemp("untemplated"), untemplated=0.3)
    return data, run_inspect(data, tmp_path_factory.mktemp("untemplated_run") / "lib")


def footprint_reads(data):
    return [r for r in data["reads"] if r[4] is not None]


def test_read_counts(clean):
    data, out = clean
    s = out["summary"]
    assert s["reads_forward"] == len(data["reads"])
    assert s["alignments_dropped_cigar"] == 5
    assert (s["refs_not_in_gtf"], s["alignments_dropped_ref_not_in_gtf"]) == (1, 9)
    assert s["reads"] == len(footprint_reads(data))
    nh = [len(r[3]) for r in footprint_reads(data)]
    assert (s["reads_nh1"], s["reads_nh2"]) == (nh.count(1), nh.count(2))
    assert out["lengths"]["reads"].sum() == pytest.approx(s["reads"])


def test_classes(clean):
    _, out = clean
    t = out["lengths"]
    assert t.loc[21, "class"] == "short state"
    assert t.loc[27, "class"] == "mixed"
    assert t.loc[27, "start_beside"] > 0.7
    assert list(t.loc[[28, 29], "class"]) == ["clean", "clean"]
    assert t.loc[20, "class"] != "short state"
    # 34 nt is clean too, but 28-29 has more reads
    assert out["summary"]["read_lengths"] == "28-29"
    assert out["summary"]["mixed_lengths"] == "27"
    assert out["summary"]["flag_mixed_lengths"] == "info"


def test_ends_at_starts_and_stops(clean):
    _, out = clean
    t = out["lengths"]
    assert list(t.loc[[21, 27, 28, 29, 34], "start_peak"]) == [-12, -12, -12, -12, -13]
    assert list(t.loc[[21, 28, 29, 34], "offset"]) == [12, 12, 12, 13]
    # 34 nt: the 5' end 1 nt upstream, so the extra length (6 nt) is 5 nt at the 3' end; 21 nt is short at the 3' end
    assert t.loc[34, "stop_peak"] - t.loc[28, "stop_peak"] == 5
    assert t.loc[28, "stop_peak"] - t.loc[21, "stop_peak"] == 7
    a = out["anchors"]
    row = a[(a["length"] == 28) & (a["anchor"] == "start")].set_index("x")["reads"]
    assert row.idxmax() == -12 and row.sum() == pytest.approx(t.loc[28, "start_reads"])


def test_fragments_by_region(clean):
    _, out = clean
    t = out["lengths"]
    assert t.loc[20, "CDS"] < 0.05 and t.loc[20, "trailer"] + t.loc[20, "transcript"] > 0.95
    assert (t.loc[[21, 27, 28, 29, 34], "CDS"] > 0.95).all()
    assert out["summary"]["length_distance"] > 0.9
    assert np.allclose(t[["leader", "CDS", "trailer", "transcript"]].sum(axis=1), 1)


def test_complexity(clean):
    data, out = clean
    refs = list(data["txs"]) + ["unannotated"]
    keys = {min((refs.index(u), p, r[4][2]) for u, p in r[3]) for r in footprint_reads(data)}
    c = out["complexity"]
    assert c.iloc[-1]["reads"] == len(footprint_reads(data))
    assert c.iloc[-1]["distinct"] == len(keys)
    assert (np.diff(c["reads"]) > 0).all() and (np.diff(c["distinct"]) > 0).all()
    assert 0.4 < c.loc[c["fraction"] == 0.5, "reads"].iloc[0] / c.iloc[-1]["reads"] < 0.6
    s = out["summary"]
    assert s["distinct_share"] == pytest.approx(len(keys) / len(footprint_reads(data)))
    assert np.isnan(s["distinct_share_at_depth"])   # fewer than COMPLEXITY_DEPTH reads


def test_expected_distinct():
    # reads per footprint 3, 1, 1; all 10 pairs of the 5 reads: 3 x 1 distinct (both from the first), 7 x 2
    c = np.array([3, 1, 1])
    assert inspect.expected_distinct(c, 2) == pytest.approx(1.7)
    assert inspect.expected_distinct(c, 1) == pytest.approx(1)
    assert inspect.expected_distinct(c, 4) == pytest.approx(2.6)   # 4 of 5: B or C left out in 2 of 5
    assert inspect.expected_distinct(c, 5) == 3
    assert np.isnan(inspect.expected_distinct(c, 6))


def test_expected_distinct_matches_the_subsamples(clean):
    data, _ = clean
    anno = annotation.load_annotation(str(data["gtf"]), str(data["fasta"]))
    curve, counts = inspect.complexity_table(bam.read_bam(str(data["bam"]), anno, {}))
    assert counts.sum() == curve["reads"].iloc[-1] and counts.size == curve["distinct"].iloc[-1]
    for n, d in zip(curve["reads"], curve["distinct"]):
        assert inspect.expected_distinct(counts, int(n)) / n == pytest.approx(d / n, abs=0.01)


def test_untemplated(clean, untemplated):
    data, out = untemplated
    reads = footprint_reads(data)
    added = [r for r in reads if len(r[1]) > r[4][2]]
    s = out["summary"]
    # every added nt mismatches the transcript, and the synthetic reads have no other mismatches (a clipped nt before
    # transcript nt 0 is not compared)
    assert s["mismatch5"] == pytest.approx(len(added) / len(reads), abs=2e-3)
    assert s["mismatch5_background"] == 0
    assert s["clip5_share"] == pytest.approx(sum(r[2].startswith("1S") for r in reads) / len(reads))
    assert s["mismatch5_T"] > 0.6
    assert (s["flag_untemplated"], s["repair"]) == ("fail", "rule 1")
    # rule 1 gives back the footprints of the run without the added nt
    t, t0 = out["lengths"], clean[1]["lengths"]
    rep = t.loc[t["reads_repaired"] > 0, ["reads_repaired", "class_repaired"]]
    assert np.allclose(rep["reads_repaired"], t0["reads"])
    assert list(rep["class_repaired"]) == list(t0["class"])
    assert s["clean_share_repaired"] == pytest.approx(clean[1]["summary"]["clean_share"])
    assert s["clean_share"] < s["clean_share_repaired"] - 0.1
    assert t.loc[29, "top_frame_repaired"] > t.loc[29, "top_frame"]


def test_clean_flags(clean):
    s = clean[1]["summary"]
    assert (s["flag_untemplated"], s["flag_end_spike"], s["repair"]) == ("pass", "pass", "no")
    assert s["flag_weak_frame"] == "pass"
    assert s["flag_short_state"] == "unknown"   # fewer than SHORT_READS reads at 20-22 nt
    assert s["flag_low_complexity"] == "unknown"   # fewer than COMPLEXITY_DEPTH reads


def test_rerun_is_identical(clean, tmp_path):
    data, out = clean
    again = run_inspect(data, tmp_path / "lib")
    for name in TABLES:
        assert filecmp.cmp(f"{out['prefix']}.{name}.tsv", f"{again['prefix']}.{name}.tsv", shallow=False), name


def test_read_bam_with_ends_same_alignments(clean):
    data, _ = clean
    anno = annotation.load_annotation(str(data["gtf"]), str(data["fasta"]))
    a = bam.read_bam(str(data["bam"]), anno, {})
    b, ends = bam.read_bam(str(data["bam"]), anno, {}, with_ends=True)
    for f in ("read", "tx", "pos5", "length", "ref"):
        assert np.array_equal(getattr(a, f), getattr(b, f))
    assert ends.clip5.size == a.read.size and ends.head.shape == (a.n_reads, bam.K)
    r = data["reads"][0]
    assert bytes(ends.head[0]).decode() == r[1][:bam.K] and bytes(ends.tail[0]).decode() == r[1][-bam.K:]


BASE = dict(classed_lengths=10, clean_share=0.6, long_share=0.02, mismatch5=0.01, mismatch5_background=0.005,
            clip5_share=0.0, mismatch5_inner=0.006, mismatch3=0.01, mismatch3_background=0.006, mixed_lengths="",
            short_reads=1e6, short_cds_share=0.8, distinct_share_at_depth=0.5)


@pytest.mark.parametrize("change, flag, state", [
    ({}, "weak_frame", "pass"),
    ({"clean_share": 0.4}, "weak_frame", "warn"),
    ({"classed_lengths": 0}, "weak_frame", "unknown"),
    ({"clean_share": 0.1}, "under_digestion", "pass"),
    ({"clean_share": 0.1, "long_share": 0.3}, "under_digestion", "fail"),
    ({"long_share": 0.3}, "under_digestion", "pass"),
    ({"mismatch5": 0.2}, "untemplated", "fail"),
    ({"mismatch5": 0.2, "clip5_share": 0.19}, "untemplated", "info"),
    ({"mismatch5": 0.03}, "untemplated", "pass"),
    ({"mismatch5": np.nan}, "untemplated", "unknown"),
    ({"mismatch5_inner": 0.05}, "end_spike", "warn"),
    ({"mismatch3": 0.05}, "end_spike", "warn"),
    ({"mixed_lengths": "27,30"}, "mixed_lengths", "info"),
    ({"short_reads": 1000}, "short_state", "unknown"),
    ({"short_cds_share": 0.3}, "short_state", "warn"),
    ({}, "short_state", "info"),
    ({}, "low_complexity", "pass"),
    ({"distinct_share_at_depth": 0.3}, "low_complexity", "warn"),
    ({"distinct_share_at_depth": 0.01}, "low_complexity", "warn"),   # no fail
    ({"distinct_share_at_depth": np.nan}, "low_complexity", "unknown"),
])
def test_flag_rules(change, flag, state):
    assert inspect.flags({**BASE, **change})[flag] == state


def test_suggest_takes_the_run_with_most_reads():
    t = pd.DataFrame({"length": [25, 26, 27, 28, 29, 31], "reads": [5, 50, 1, 60, 70, 500],
                      "class": ["clean", "clean", "mixed", "clean", "clean", "clean"]})
    assert inspect.suggest(t) == "31-31"
    t.loc[5, "reads"] = 100
    assert inspect.suggest(t) == "28-29"
    assert inspect.suggest(t.assign(**{"class": "no frame"})) == ""

import filecmp

import numpy as np
import pandas as pd
import pysam
import pytest

from pitstop import inspect, repair

from synth import make_inspect_dataset


def run_repair(data, prefix, bam=None):
    repair.run(str(bam or data["bam"]), str(data["gtf"]), str(data["fasta"]), str(prefix))
    return prefix


def records(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as f:
        return list(f.fetch(until_eof=True))


def header(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as f:
        return f.header.to_dict()


@pytest.fixture(scope="module")
def untemplated(tmp_path_factory):
    d = tmp_path_factory.mktemp("untemplated")
    data = make_inspect_dataset(d, untemplated=0.3)
    return data, run_repair(data, d / "lib.repaired")


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
    d = tmp_path_factory.mktemp("clean")
    data = make_inspect_dataset(d)
    return data, run_repair(data, d / "lib.repaired")


def test_untemplated_nt_are_clipped(untemplated):
    data, prefix = untemplated
    aligned = {name: foot for name, _, _, _, foot, al in data["reads"] if al}   # the added nt is in the alignment
    before, after = records(data["bam"]), records(f"{prefix}.bam")
    assert len(before) == len(after)
    n = 0
    for a, b in zip(before, after):
        if a.query_name in aligned:
            n += 1
            assert b.cigarstring == f"1S{aligned[a.query_name][2]}M"
            assert b.reference_start == a.reference_start + 1
            assert b.query_sequence == a.query_sequence
        else:
            assert b.to_string() == a.to_string()
    assert n > 1000
    t = pd.read_csv(f"{prefix}.repair.tsv", sep="\t")
    assert t.loc[t["clipped"] == 1, "reads"].sum() == len(aligned)
    assert set(t["clipped"]) == {0, 1}


def test_inspect_predicts_the_repair(untemplated, tmp_path):
    data, prefix = untemplated
    gtf, fasta = str(data["gtf"]), str(data["fasta"])
    s0 = inspect.run(str(data["bam"]), gtf, fasta, str(tmp_path / "before"), report=False)
    s1 = inspect.run(f"{prefix}.bam", gtf, fasta, str(tmp_path / "after"), report=False)
    t0 = pd.read_csv(tmp_path / "before.lengths.tsv", sep="\t")
    t1 = pd.read_csv(tmp_path / "after.lengths.tsv", sep="\t")
    t0, t1 = t0[t0["reads_repaired"] > 0], t1[t1["reads"] > 0]
    assert list(t0["length"]) == list(t1["length"])
    assert np.allclose(t0["reads_repaired"], t1["reads"])
    assert np.allclose(t0["top_frame_repaired"], t1["top_frame"], equal_nan=True)
    assert list(t0["class_repaired"]) == list(t1["class"])
    assert s1["clean_share"] == pytest.approx(s0["clean_share_repaired"])
    assert s0["flag_untemplated"] == "fail" and s0["repair"] == "yes"
    assert s1["flag_untemplated"] == "info" and s1["repair"] == "no"   # now a 5' soft clip


def test_second_run_changes_nothing(untemplated, tmp_path):
    data, prefix = untemplated
    again = run_repair(data, tmp_path / "again", bam=f"{prefix}.bam")
    assert [a.to_string() for a in records(f"{again}.bam")] == [a.to_string() for a in records(f"{prefix}.bam")]
    assert set(pd.read_csv(f"{again}.repair.tsv", sep="\t")["clipped"]) == {0}
    pg = header(f"{again}.bam")["PG"]
    assert [p["ID"] for p in pg] == ["pitstop", "pitstop.1"] and pg[1]["PP"] == "pitstop"


def test_clean_bam_changes_only_the_header(clean):
    data, prefix = clean
    assert [a.to_string() for a in records(f"{prefix}.bam")] == [a.to_string() for a in records(data["bam"])]
    assert set(pd.read_csv(f"{prefix}.repair.tsv", sep="\t")["clipped"]) == {0}
    h0, h1 = header(data["bam"]), header(f"{prefix}.bam")
    assert [p["PN"] for p in h1.pop("PG")] == ["pitstop"] and "PG" not in h0
    assert h0 == h1


def test_rerun_is_byte_identical(untemplated, tmp_path):
    data, prefix = untemplated
    again = run_repair(data, tmp_path / prefix.name)
    for suffix in ("bam", "repair.tsv"):
        assert filecmp.cmp(f"{prefix}.{suffix}", f"{again}.{suffix}", shallow=False), suffix


def test_sorted_bam_stays_sorted(untemplated, tmp_path):
    data, prefix = untemplated
    sorted_bam = tmp_path / "sorted.bam"
    pysam.sort("--no-PG", "-o", str(sorted_bam), str(data["bam"]))
    out = run_repair(data, tmp_path / "sorted.repaired", bam=sorted_bam)
    assert header(f"{out}.bam")["HD"]["SO"] == "coordinate"
    after = records(f"{out}.bam")
    keys = [(a.reference_id if a.reference_id >= 0 else 1 << 30, a.reference_start) for a in after]
    assert keys == sorted(keys)
    # the same alignments as from the unsorted BAM, and the clips moved some past others
    unsorted = {a.to_string() for a in records(f"{prefix}.bam")}
    assert sorted(a.to_string() for a in after) == sorted(unsorted)
    moved = {(a.query_name, a.flag, a.reference_id, a.reference_start): b.reference_start
             for a, b in zip(records(data["bam"]), records(f"{prefix}.bam"))}
    naive = [(a.reference_id, moved[a.query_name, a.flag, a.reference_id, a.reference_start])
             for a in records(sorted_bam)]
    assert naive != sorted(naive)
    with pysam.AlignmentFile(f"{out}.bam", "rb") as f:    # the index
        ref = f.references[0]
        assert sum(1 for _ in f.fetch(ref)) == sum(1 for a in after if a.reference_name == ref)


@pytest.mark.parametrize("md, k, out, mism", [
    ("0T29", 1, "29", 1), ("0T0C29", 2, "29", 2), ("0T0C0A29", 2, "0A29", 2), ("30", 1, "29", 0),
    ("1T27", 2, "27", 1), ("0T3^AC20", 1, "3^AC20", 1)])
def test_clip_md(md, k, out, mism):
    assert repair.clip_md(md, k) == (out, mism)


def test_rule_one_cases(tmp_path):
    """Hand-made alignments on one transcript: the soft clip, the position, NM and MD, and what repair leaves."""
    data = make_inspect_dataset(tmp_path)
    txs = data["txs"]
    t = next(n for n, tx in txs.items() if tx["cds_start"] is not None)
    S, p, L = txs[t]["seq"], 50, 29
    f = S[p:p + L]

    def x(i):   # a nt that does not match S[i]
        return "C" if S[i] == "A" else "A"

    nm_md = lambda nm, md: [("NH", 1, "C"), ("NM", nm, "C"), ("MD", md, "Z"), ("XA", 0, "i")]
    cases = [  # name, flag, ref, pos, CIGAR, sequence, tags -> CIGAR, pos, NM, MD after the repair
        ("one", 0, t, p - 1, "30M", x(p - 1) + f, nm_md(1, f"0{S[p - 1]}29"), "1S29M", p, 0, "29"),
        ("two", 0, t, p - 2, "31M", x(p - 2) + x(p - 1) + f, nm_md(2, f"0{S[p - 2]}0{S[p - 1]}29"),
         "2S29M", p, 0, "29"),
        ("three", 0, t, p - 3, "32M", x(p - 3) + x(p - 2) + x(p - 1) + f,
         nm_md(3, f"0{S[p - 3]}0{S[p - 2]}0{S[p - 1]}29"), "2S30M", p - 1, 1, f"0{S[p - 1]}29"),
        ("clipped1", 0, t, p - 1, "1S30M", "G" + x(p - 1) + f, nm_md(1, f"0{S[p - 1]}29"), "2S29M", p, 0, "29"),
        ("clipped2", 0, t, p - 1, "2S30M", "GG" + x(p - 1) + f, nm_md(1, f"0{S[p - 1]}29"), None, None, None, None),
        ("eqx", 0, t, p - 1, "1X29=", x(p - 1) + f, [], "1S29=", p, None, None),
        ("reverse", 16, t, p - 1, "30M", x(p - 1) + f, [], None, None, None, None),
        ("multi", 0, t, p - 1, "30M", x(p - 1) + f, [("NM", 1, "C")], "1S29M", p, 0, None),
        ("multi", 256, t, p - 1, "30M", None, [("NM", 1, "C")], "1S29M", p, 0, None),   # no SEQ: the read's before
        ("nseq", 0, t, p - 1, "30M", "N" + f, [], None, None, None, None),
        ("unknown", 0, "unannotated", 5, "30M", "C" * 30, [], None, None, None, None),
        ("lonely", 256, t, p - 1, "30M", None, [], None, None, None, None),       # no SEQ, and not after its read
        ("short", 0, t, p - 1, "2M", x(p - 1) + x(p), [], "1S1M", p, None, None),  # one nt stays aligned
    ]
    refs = list(txs) + ["unannotated"]
    head = {"HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": r, "LN": len(txs[r]["seq"]) if r in txs else 1000} for r in refs]}
    path = tmp_path / "cases.bam"
    with pysam.AlignmentFile(str(path), "wb", header=head) as out:
        for name, flag, ref, pos, cigar, seq, tags, *_ in cases:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, refs.index(ref), pos
            a.mapping_quality, a.cigarstring, a.query_sequence = 255, cigar, seq
            a.set_tags(tags)
            out.write(a)
    prefix = run_repair(data, tmp_path / "cases.repaired", bam=path)
    for case, a, b in zip(cases, records(path), records(f"{prefix}.bam")):
        name, cigar, pos, nm, md = case[0], *case[7:]
        if cigar is None:
            assert b.to_string() == a.to_string(), name
            continue
        assert (b.cigarstring, b.reference_start) == (cigar, pos), name
        assert b.query_sequence == a.query_sequence, name
        assert [g[0] for g in b.get_tags()] == [g[0] for g in a.get_tags()], name   # the order of the tags
        if nm is not None:
            assert b.get_tag("NM") == nm, name
        if md is not None:
            assert b.get_tag("MD") == md, name
    t = pd.read_csv(f"{prefix}.repair.tsv", sep="\t")
    assert t.values.tolist() == [[2, 1, 1, 1], [30, 0, 2, 2], [30, 1, 5, 4], [31, 2, 1, 1], [32, 2, 1, 1]]

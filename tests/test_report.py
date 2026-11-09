import filecmp
import html
import json
import pathlib
import re

import numpy as np
import pandas as pd
import pytest

from pitstop import inspect, report, starts

from synth import make_inspect_dataset, make_start_dataset
from test_starts import pool_sheet

MQC = ("inspect_general_mqc.json", "inspect_lengths_mqc.json", "inspect_frame_mqc.json")


def run_inspect(data, prefix):
    inspect.run(str(data["bam"]), str(data["gtf"]), str(data["fasta"]), str(prefix))
    return prefix


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    out = {}
    for name, u in (("clean", 0.0), ("untemplated", 0.3)):
        data = make_inspect_dataset(tmp_path_factory.mktemp(name), untemplated=u)
        out[name] = (data, run_inspect(data, tmp_path_factory.mktemp(f"{name}_run") / "lib"))
    return out


def page(prefix):
    with open(f"{prefix}.inspect.html", encoding="utf-8") as f:
        return f.read()


def visible_text(doc):
    """The text a reader sees outside the figures (and outside the code of the page)."""
    doc = re.sub(r"<(svg|style|script)\b.*?</\1>", " ", doc, flags=re.S)
    return html.unescape(re.sub(r"<[^>]+>", " ", doc))


def allowed_numbers(prefix, tables=("summary", "lengths", "anchors", "ends", "complexity"),
                    modules=(inspect, report.inspect, report.base)):
    """Every number of the tables, and every constant of the modules, as the page may format it."""
    values, strings = [], []
    for name in tables:
        if not pathlib.Path(f"{prefix}.{name}.tsv").exists():
            continue
        df = pd.read_csv(f"{prefix}.{name}.tsv", sep="\t", dtype=str, keep_default_na=False)
        for v in df.to_numpy().ravel():
            try:
                values.append(float(v))
            except ValueError:
                strings.append(v)
    for mod in modules:
        for k, v in vars(mod).items():
            if k.isupper() and isinstance(v, (int, float, tuple)) and not isinstance(v, bool):
                values += [x for x in (v if isinstance(v, tuple) else (v,)) if isinstance(x, (int, float))]
    out = set()
    for v in values:
        if np.isnan(v):
            continue
        for x in (v, 100 * v):
            out.update(f"{x:.{d}f}" for d in range(5))
            out.add(f"{x:,.0f}")
            out.add(str(x))
    out.update(str(i) for i in range(-100, 101))
    return out, " ".join(strings)


@pytest.fixture(scope="module")
def under_digested(runs, tmp_path_factory):
    """The clean run's tables with a 3' end spread of a failing library (the synthetic data have too few reads at the
    starts to measure it), and the report drawn from them."""
    _, src = runs["clean"]
    prefix = tmp_path_factory.mktemp("under_digested") / "lib"
    for name in ("summary", "lengths", "anchors", "ends", "complexity"):
        df = pd.read_csv(f"{src}.{name}.tsv", sep="\t", dtype=str, keep_default_na=False)
        if name == "summary":
            df[["end3_reads", "end3_spread", "flag_under_digestion"]] = ["5000", "1.9123", "fail"]
        df.to_csv(f"{prefix}.{name}.tsv", sep="\t", index=False)
    report.inspect_report(str(prefix))
    return prefix


def test_under_digestion(runs, under_digested):
    text = visible_text(page(under_digested))
    assert "repair does not help with under-digestion." in text and "spread 1.91 nt" in text
    assert "under-digestion. No rule per read fixes this." in text
    clean = visible_text(page(runs["clean"][1]))
    assert "too few reads at the starts" in clean and "under-digestion." not in clean


@pytest.mark.parametrize("name", ["clean", "untemplated", "under_digested"])
def test_every_number_is_in_a_table(runs, under_digested, name):
    prefix = under_digested if name == "under_digested" else runs[name][1]
    allowed, strings = allowed_numbers(prefix)
    text = visible_text(page(prefix)).replace("−", "-")
    tokens = re.findall(r"(?<![\w.\-])-?\d[\d,]*(?:\.\d+)?(?![\w])", text)
    missing = {t for t in tokens if t.rstrip(",") not in allowed and t.rstrip(",") not in strings}
    assert not missing


def test_rerun_is_byte_identical(runs, tmp_path):
    data, prefix = runs["untemplated"]
    again = run_inspect(data, tmp_path / "lib")
    for suffix in ("inspect.html",) + MQC:
        assert filecmp.cmp(f"{prefix}.{suffix}", f"{again}.{suffix}", shallow=False), suffix


def test_no_paths_and_no_network(runs):
    data, prefix = runs["untemplated"]
    doc = page(prefix)
    for path in (prefix.parent, data["bam"].parent):
        assert str(path) not in doc
    assert not re.search(r"src=\"https?:", doc)   # loads nothing
    links = re.findall(r"href=\"(https?:[^\"]*)\"", doc)
    assert links and all(u.startswith(report.DOCS) for u in links)   # links only to the docs of the flags
    assert "reads.bam" in doc and "annotation.gtf" in doc   # file names only


def test_flag_pages_exist(runs):
    """Each flag of the report links to its page, and the docs have that page."""
    doc = page(runs["clean"][1])
    pages = re.findall(rf"href=\"{re.escape(report.DOCS)}([^\"/]+)/\"", doc)
    assert sorted(pages) == sorted(f.replace("_", "-") for f in inspect.FLAGS)
    docs = pathlib.Path(__file__).resolve().parent.parent / "docs" / "flags"
    for p in pages:
        assert (docs / f"{p}.md").is_file(), p


def test_flags_and_steps(runs):
    _, prefix = runs["untemplated"]
    text = visible_text(page(prefix))
    for name in ("Weak frame", "Long reads", "Under-digestion", "Untemplated 5′ nt", "End spike", "Mixed lengths",
                 "Short state", "Low complexity"):
        assert name in text
    assert "For the next alignment: let the aligner clip the added nucleotide." in text and "bowtie2 --local" in text
    assert text.index("For the next alignment") < text.index("For this BAM: run repair.")
    clean = visible_text(page(runs["clean"][1]))
    assert "Use the clean lengths: 28–29 nt." in clean and "run repair" not in clean and "For the next" not in clean
    assert "Mixed lengths: 27, 29 nt." in clean and "Their frame in the CDS decides" in clean


def test_figures_use_theme_colours(runs):
    doc = page(runs["clean"][1])
    figures = re.findall(r"<svg role=\"img\".*?</svg>", doc, flags=re.S)
    assert len(figures) == 8
    for fig in figures:
        assert not re.search(r"#[0-9a-fA-F]{6}\b", fig)
        assert "var(--" in fig
    ids = re.findall(r'\bid="([^"]+)"', doc)
    assert len(ids) == len(set(ids))


def test_multiqc_files(runs):
    _, prefix = runs["clean"]
    general, lengths, frame = (json.load(open(f"{prefix}.{m}")) for m in MQC)
    assert general["plot_type"] == "generalstats" and list(general["data"]) == ["lib"]
    row = general["data"]["lib"]
    s = pd.read_csv(f"{prefix}.summary.tsv", sep="\t").iloc[0]
    assert row["clean_share"] == pytest.approx(s["clean_share"], abs=1e-4)
    assert lengths["plot_type"] == frame["plot_type"] == "linegraph"
    t = pd.read_csv(f"{prefix}.lengths.tsv", sep="\t")
    assert [x for x, _ in lengths["data"]["lib"]] == t["length"].tolist()
    assert all(0 <= y <= 1 for _, y in frame["data"]["lib"])


def test_complexity_at_depth(runs, tmp_path, monkeypatch):
    data, _ = runs["clean"]
    monkeypatch.setattr(inspect, "COMPLEXITY_DEPTH", 1000)
    prefix = run_inspect(data, tmp_path / "lib")
    s = pd.read_csv(f"{prefix}.summary.tsv", sep="\t").iloc[0]
    assert s["flag_low_complexity"] == "pass" and inspect.COMPLEXITY < s["distinct_share_at_depth"] < 1
    text = visible_text(page(prefix))
    assert f"In a subsample of 1,000 reads, {s['distinct_share_at_depth']:.2f} are distinct footprints." in text
    allowed, strings = allowed_numbers(prefix)
    tokens = re.findall(r"(?<![\w.\-])-?\d[\d,]*(?:\.\d+)?(?![\w])", text.replace("−", "-"))
    assert not {t for t in tokens if t.rstrip(",") not in allowed and t.rstrip(",") not in strings}
    general = json.load(open(f"{prefix}.inspect_general_mqc.json"))
    assert general["data"]["lib"]["distinct_share_at_depth"] == pytest.approx(s["distinct_share_at_depth"], abs=1e-4)


# ---------- the report of starts ----------

START_TABLES = ("starts_stats", "lengths", "libraries", "pairs", "start_codons", "false_calls")


@pytest.fixture(scope="module")
def started(tmp_path_factory):
    """starts on synthetic pairs: with an ORF table (s), a scan (scan), and a window without a length (fail: no start
    lengths, the run stops after its report)."""
    d = tmp_path_factory.mktemp("startreport")
    data = make_start_dataset(d)
    sheet = pool_sheet(d, data)
    args = (str(sheet), str(data["gtf"]), str(data["fasta"]))
    starts.run(*args, str(d / "s"), str(data["orfs"]))
    starts.run(*args, str(d / "scan"), scan=True)
    with pytest.raises(ValueError, match="no read length passes the QC"):
        starts.run(*args, str(d / "fail"), str(data["orfs"]), window=(31, 34))
    return dict(dir=d, data=data, sheet=sheet, args=args)


def starts_page(prefix):
    with open(f"{prefix}.starts.html", encoding="utf-8") as f:
        return f.read()


@pytest.mark.parametrize("name", ["s", "scan", "fail"])
def test_starts_every_number_is_in_a_table(started, name):
    prefix = started["dir"] / name
    allowed, strings = allowed_numbers(prefix, START_TABLES, (starts, report.starts, report.base))
    text = visible_text(starts_page(prefix)).replace("−", "-")
    tokens = re.findall(r"(?<![\w.\-])-?\d[\d,]*(?:\.\d+)?(?![\w])", text)
    missing = {t for t in tokens if t.rstrip(",") not in allowed and t.rstrip(",") not in strings}
    assert not missing


def test_starts_rerun_is_byte_identical(started, tmp_path):
    d = started["dir"]
    starts.run(*started["args"], str(tmp_path / "scan"), scan=True)
    assert filecmp.cmp(d / "scan.starts.html", tmp_path / "scan.starts.html", shallow=False)


def test_starts_no_paths_and_no_network(started):
    for name in ("s", "scan", "fail"):
        doc = starts_page(started["dir"] / name)
        assert str(started["dir"]) not in doc and str(started["data"]["gtf"].parent) not in doc
        assert not re.search(r"(src|href)=\"https?:", doc)          # loads nothing, links nowhere
        assert "harringtonine1" in doc and "annotation.gtf" in doc     # file names only


def test_starts_report_sections(started):
    text = visible_text(starts_page(started["dir"] / "scan"))
    titles = re.findall(r'<h2 id="t-[^"]+">([^<]+)</h2>', starts_page(started["dir"] / "scan"))
    assert titles == ["Verdict and next steps", "Calls", "Pairs", "Start lengths", "Pair detail", "Not in this report"]
    assert "Per pool: harringtonine at the annotated starts" in text
    assert "Start lengths: 28–30 nt." in text and "Every pair has a read-length distance below 0.4" in text
    assert "In every pair the initiation library has the higher start peak." in text
    assert "5′ UTR: ATG is" in text and "catalogue.tsv has the uORFs" in text
    assert "pitstop starts --libraries sheet.tsv --gtf annotation.gtf --fasta genome.fa --scan --out-prefix scan" in text


def test_starts_pair_picker(started):
    doc = starts_page(started["dir"] / "scan")
    options = re.findall(r'<option value="(pair-\d+)"', doc)
    panels = re.findall(r'<div class="libpanel" id="(pair-\d+)"[^>]*?( hidden)?>', doc)
    assert len(options) == len(panels) == 2 and options == [p for p, _ in panels]
    assert [h for _, h in panels] == ["", " hidden"]                  # the first pair shows
    assert doc.index('id="pair-pick"') < doc.index('id="pair-0"')    # the picker is above the panels it controls
    assert re.findall(r'data-pair="(pair-\d+)"', doc) == options     # the side bar picks the same pairs
    figures = re.findall(r"<svg role=\"img\".*?</svg>", doc, flags=re.S)
    for fig in figures:
        assert not re.search(r"#[0-9a-fA-F]{6}\b", fig) and "var(--" in fig
    ids = re.findall(r'\bid="([^"]+)"', doc)
    assert len(ids) == len(set(ids))


def test_starts_report_when_qc_fails(started):
    text = visible_text(starts_page(started["dir"] / "fail"))
    assert "No read length passes the QC, so starts stopped." in text and "--start-lengths LO-HI" in text
    assert "not computed, no start lengths." in text and "No calls: no start lengths." in text
    assert "Read-length distance" in text                              # the pairs do not need start lengths


@pytest.fixture(scope="module")
def bad_pairs(started, tmp_path_factory):
    """The scan's tables with pair 2 of another digestion and pair 1 without a harringtonine effect, and the report
    drawn from them."""
    src = started["dir"] / "scan"
    prefix = tmp_path_factory.mktemp("bad_pairs") / "scan"
    for name in ("starts_stats", "lengths", "profiles", "libraries", "pairs", "kernel", "start_codons",
                 "false_calls"):
        df = pd.read_csv(f"{src}.{name}.tsv", sep="\t", dtype=str, keep_default_na=False)
        if name == "pairs":
            df.loc[1, "distance"] = "0.6123"
        if name == "libraries":
            harr = df["library"].str.endswith("harringtonine1")
            df.loc[harr, "start_enrichment"] = "1.5"
        df.to_csv(f"{prefix}.{name}.tsv", sep="\t", index=False)
    report.starts_report(str(prefix))
    return prefix


def test_starts_bad_pairs(bad_pairs):
    text = visible_text(starts_page(bad_pairs))
    assert "Check the pairs of another digestion: A · pair 2." in text and "Read-length distance 0.61:" in text
    assert "Harringtonine did not raise the start peak: A · pair 1." in text
    assert "harringtonine did not raise it. Leave the pair out, or make the library again." in text
    allowed, strings = allowed_numbers(bad_pairs, START_TABLES, (starts, report.starts, report.base))
    tokens = re.findall(r"(?<![\w.\-])-?\d[\d,]*(?:\.\d+)?(?![\w])", text.replace("−", "-"))
    assert not {t for t in tokens if t.rstrip(",") not in allowed and t.rstrip(",") not in strings}


def test_distance_near_the_cut():
    # a distance that 2 decimals would round onto the cut shows 4, so a pass never reads as the cut
    assert report.starts.dist(0.39991) == "0.3999" and report.starts.dist(0.41) == "0.41"
    assert report.starts.dist(starts.PAIR_DISTANCE) == "0.40" and report.starts.dist(0.1) == "0.10"

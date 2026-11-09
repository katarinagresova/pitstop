"""The report of `inspect` (<prefix>.inspect.html), one library, from the tables of its run; and its MultiQC files."""
import json
import logging

import numpy as np
import pandas as pd

from .. import inspect as qc
from ..logo import LOCKUP
from .base import (C, HEAT_BINS, RC, command, count, esc, f2, fig_block, findings, html_table, icon, int_list, isna,
                   legend, nt, page, pct, section, span, state_cell, svg, tiles, tsv_view)

log = logging.getLogger(__name__)

WIDE = (8.0, 2.6)          # inches; the page scales a figure to its column
SMALL = (3.9, 2.3)
# The cuts of the lines of text that are not flags (each line says its number)
CDS_LOW = 0.7        # a lower CDS share gets the advice on contamination
DISTANCE = 0.3       # a larger length distance: the reads without a CDS are mostly not ribosomes
RISING = 0.3         # new footprints per read at full depth at this or more: the curve still rises
BIAS = 1.5           # an end preference this strong (or 1 / BIAS) gets advice
SHOWN = 0.001        # the figures show the lengths from the first to the last with this share of the reads
CHANGED = 0.02       # the frame section names the clean share after repair if it changes by this much
HIDDEN = 0.005       # the next steps name the hidden added nt from this share
DOCS = "https://katarinagresova.github.io/pitstop/flags/"   # a flag's page: DOCS + its name with "-" for "_"


TEXT = (["library", "pitstop", "bam", "gtf", "fasta", "clean_lengths", "read_lengths", "clean_lengths_repaired",
         "read_lengths_repaired", "mixed_lengths", "bias_max_at", "repair"] + [f"flag_{f}" for f in qc.FLAGS])
CLASSES = (("clean", "clean"), ("short state", "short state"), ("no frame", "no frame"), ("few reads", "too few reads"))
SECTIONS = (("verdict", "Verdict and next steps"), ("reads", "The reads that pitstop reads"),
            ("depth", "Depth and complexity"), ("lengths", "Read lengths"), ("regions", "Regions per length"),
            ("frame", "Frame per length"), ("ends", "Ends at the annotated starts and stops"),
            ("mismatch", "Read ends against the transcript"), ("scope", "Not in this report"))


def load(prefix):
    """The tables of an inspect run: summary (a dict: text as str, numbers as float), lengths, anchors, ends,
    complexity."""
    raw = pd.read_csv(f"{prefix}.summary.tsv", sep="\t", dtype=str, keep_default_na=False).iloc[0].to_dict()
    s = {k: (v if k in TEXT else float(v) if v not in ("", "NA") else np.nan) for k, v in raw.items()}
    lengths = pd.read_csv(f"{prefix}.lengths.tsv", sep="\t")
    anchors = pd.read_csv(f"{prefix}.anchors.tsv", sep="\t")
    ends = pd.read_csv(f"{prefix}.ends.tsv", sep="\t", dtype={"length": str})
    cplx = pd.read_csv(f"{prefix}.complexity.tsv", sep="\t")
    return s, lengths, anchors, ends, cplx


def shown_lengths(t):
    """The lengths the figures draw: from the first to the last with 0.1% of the reads or more."""
    keep = t.loc[t["share"] >= SHOWN, "length"]
    if keep.empty:
        keep = t["length"]
    return t[(t["length"] >= keep.min()) & (t["length"] <= keep.max())].reset_index(drop=True)


def bar_axis(ax, L, ylim=None, ylabel=None, pct_axis=True):
    ax.set_xlim(L.min() - 0.7, L.max() + 0.7)
    step = 1 if len(L) <= 26 else 2
    ax.set_xticks([x for x in L if x % step == 0])
    ax.set_xlabel("footprint length (nt)")
    if ylim:
        ax.set_ylim(*ylim)
    if ylabel:
        ax.set_ylabel(ylabel)
    if pct_axis:
        ax.yaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%")
    ax.grid(axis="y")
    ax.set_axisbelow(True)


def lengths_chart(t, lib):
    from matplotlib.figure import Figure
    fig = Figure(figsize=WIDE)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.17)
    ax = fig.add_subplot()
    L = t["length"].to_numpy()
    tips = {}
    for i, r in t.iterrows():
        cls = r["class"]
        kw = dict(width=0.78, linewidth=0.8, edgecolor=C["surface"])
        if cls == "clean":
            kw["color"] = C["k-clean"]
        elif cls == "short state":
            kw["color"] = C["k-short"]
        elif cls == "no frame":
            kw["color"] = C["k-noframe"]
        else:
            kw.update(facecolor="none", edgecolor=C["muted"], linestyle=(0, (2, 1.5)), linewidth=0.7)
        b = ax.bar([r["length"]], [r["share"]], **kw)
        gid = f"len-{i}"
        b[0].set_gid(gid)
        tips[gid] = (f"{r['length']} nt: {pct(r['share'])} of the reads ({count(r['reads'])}), {cls}"
                     + (", mixed at the starts" if r["mixed"] else ""))
    bar_axis(ax, L, (0, max(t["share"].max() * 1.12, 1e-6)), "share of the reads")
    return svg(fig, "lengths", f"Reads per footprint length of {lib}, by class", tips)


def regions_chart(t, lib):
    from matplotlib.figure import Figure
    fig = Figure(figsize=WIDE)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.17)
    ax = fig.add_subplot()
    colours = [C["r-5utr"], C["r-cds"], C["r-3utr"], C["r-nocds"]]
    labels = ["5′ UTR", "CDS", "3′ UTR", "no CDS"]
    tips = {}
    for i, r in t.iterrows():
        bottom = 0.0
        for j, (reg, c) in enumerate(zip(qc.REGIONS, colours)):
            v = 0.0 if isna(r[reg]) else r[reg]
            if v <= 0:
                continue
            b = ax.bar([r["length"]], [v], bottom=bottom, width=0.78, color=c, edgecolor=C["surface"], linewidth=0.8)
            gid = f"reg-{i}-{j}"
            b[0].set_gid(gid)
            tips[gid] = f"{r['length']} nt, {labels[j]}: {pct(v)} of its P-sites"
            bottom += v
    bar_axis(ax, t["length"].to_numpy(), (0, 1), "share of the P-sites")
    return svg(fig, "regions", f"P-sites by region per footprint length of {lib}", tips)


def frame_chart(t, lib):
    from matplotlib.figure import Figure
    fig = Figure(figsize=WIDE)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.17)
    ax = fig.add_subplot()
    tips = {}
    for i, r in t.iterrows():
        if r["interior"] <= 0 or isna(r["top_phase"]):
            continue
        top = int(r["top_phase"])
        bottom = 0.0
        for p in [top] + [p for p in range(3) if p != top]:
            v = r[f"phase{p}"]
            b = ax.bar([r["length"]], [v], bottom=bottom, width=0.78, color=C[f"f{p}"], edgecolor=C["surface"],
                       linewidth=0.8, alpha=1 if r["class"] != "few reads" else 0.35)
            gid = f"frm-{i}-{p}"
            b[0].set_gid(gid)
            tips[gid] = f"{r['length']} nt, phase {p}: {pct(v)} of {count(r['interior'])} reads in the CDS interior"
            bottom += v
        rep = r["top_frame_repaired"]
        if not isna(rep) and abs(rep - r["top_frame"]) >= 0.01:
            ax.scatter([r["length"]], [rep], s=22, facecolor=C["surface"], edgecolor=C["ink"], linewidth=1.2, zorder=4)
    ax.axhline(qc.CLEAN, color=C["ink2"], linewidth=1, linestyle=(0, (4, 3)), zorder=3)
    bar_axis(ax, t["length"].to_numpy(), (0, 1), "share of the CDS-interior reads")
    return svg(fig, "frame", f"5' phase of the CDS-interior reads per footprint length of {lib}", tips)


def depth_chart(c, lib):
    from matplotlib.figure import Figure
    fig = Figure(figsize=(8.0, 2.8))
    fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.16)
    ax = fig.add_subplot()
    n, d = c["reads"].to_numpy(dtype=float), c["distinct"].to_numpy(dtype=float)
    top = n.max()
    scale, unit = (1e6, "M") if top >= 1e6 else (1e3, "k") if top >= 1e3 else (1, "")
    ax.plot([0, top / scale], [0, top / scale], color=C["ink2"], linewidth=1, linestyle=(0, (4, 3)))
    ax.plot(np.r_[0, n] / scale, np.r_[0, d] / scale, color=C["ink"], linewidth=2)
    ax.scatter(n / scale, d / scale, s=16, color=C["ink"], edgecolor=C["surface"], linewidth=1, zorder=3)
    ax.set_xlim(0, top / scale * 1.02)
    ax.set_ylim(0, top / scale * 1.02)
    ax.set_xlabel(f"reads{' (' + unit + ')' if unit else ''}")
    ax.set_ylabel(f"distinct footprints{' (' + unit + ')' if unit else ''}")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    return svg(fig, "depth", f"Distinct footprints against the reads, {lib}")


def heat_chart(a, t, anchor, lib):
    """Rows: the shown lengths; cells: the share of the row's reads in the window (anchors.tsv), binned."""
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.figure import Figure
    sub = a[(a["anchor"] == anchor) & a["length"].isin(t["length"])]
    grid = sub.pivot(index="length", columns="x", values="share").reindex(t["length"]).fillna(0)
    xs, L = grid.columns.to_numpy(), grid.index.to_numpy()
    h = 0.6 + 0.13 * len(L)
    fig = Figure(figsize=(3.9, h))
    fig.subplots_adjust(left=0.16, right=0.98, top=1 - 0.08 / h, bottom=0.42 / h)
    ax = fig.add_subplot()
    cmap = ListedColormap([C[f"h{i}"] for i in range(8)])
    norm = BoundaryNorm((0,) + HEAT_BINS + (1.0001,), 8)
    rows = np.arange(len(L))   # one row per length, also where lengths are missing
    ax.pcolormesh(np.r_[xs - 0.5, xs[-1] + 0.5], np.r_[rows - 0.5, rows[-1] + 0.5], grid.to_numpy(), cmap=cmap,
                  norm=norm, edgecolor=C["surface"], linewidth=0.3)
    ax.set_ylim(rows[-1] + 0.5, -0.5)
    ax.set_yticks(rows, [f"{x}" + (" mixed" if m else "") for x, m in zip(L, t["mixed"])])
    ax.set_xticks([x for x in xs if x % 3 == 0], [nt(x) for x in xs if x % 3 == 0])
    ax.axvline(0, color=C["ink2"], linewidth=0.8, linestyle=(0, (2, 2)))
    ax.set_xlabel("5' end, nt from the start codon" if anchor == "start" else "3' end, nt from the stop codon")
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    return svg(fig, f"heat-{anchor}", f"{'5' if anchor == 'start' else '3'}' ends around the "
                                     f"{anchor} codons per length, {lib}")


def mismatch_chart(row, end, lib):
    from matplotlib.figure import Figure
    fig = Figure(figsize=SMALL)
    fig.subplots_adjust(left=0.17, right=0.97, top=0.93, bottom=0.2)
    ax = fig.add_subplot()
    v = np.array([row[f"mm{i}"] for i in range(qc.K)], dtype=float)
    bg = np.nanmean(v[qc.BACKGROUND[0]:qc.BACKGROUND[1]])
    tips = {}
    for i in range(qc.K):
        b = ax.bar([i], [0 if np.isnan(v[i]) else v[i]], width=0.72, color=C["ink2"], edgecolor=C["surface"],
                   linewidth=0.8)
        b[0].set_gid(f"mm{end}-{i}")
        tips[f"mm{end}-{i}"] = f"position {i}: {pct(v[i], 2)} mismatched"
    if not np.isnan(bg):
        ax.axhline(bg, color=C["muted"], linewidth=1, linestyle=(0, (4, 3)))
    top = np.nanmax(np.r_[v, 0.01])
    ax.set_ylim(0, top * 1.15)
    ax.set_xticks(range(qc.K))
    ax.set_xlabel(f"position from the {end}' end of the read")
    ax.yaxis.set_major_formatter(lambda x, _: f"{100 * x:.0f}%" if top >= 0.05 else f"{100 * x:.1f}%")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    return svg(fig, f"mm{end}", f"Mismatch rate by position from the {end}' end, {lib}", tips)


def flag_rows(s, t):
    """(id, name, state, value, cut, section) per flag."""
    st = {f: s[f"flag_{f}"] for f in qc.FLAGS}
    w, lr, ud, u, sp = qc.WEAK_FRAME, qc.LONG_READS, qc.UNDER_DIGESTION, qc.UNTEMPLATED, qc.SPIKE
    depth = count(qc.COMPLEXITY_DEPTH)
    m5, bg5 = s["mismatch5"], s["mismatch5_background"]
    ratio = m5 / bg5 if bg5 and not isna(bg5) and bg5 > 0 else np.inf
    rtext = "∞" if np.isinf(ratio) else f"{ratio:.1f}"
    if st["untemplated"] == "fail":
        uv = f"{pct(m5)}, {rtext}× background, not clipped"
    elif st["untemplated"] == "info":
        uv = f"{pct(m5)}, clipped by the aligner"
    elif st["untemplated"] == "unknown":
        uv = "no read sequences"
    else:
        uv = f"{pct(m5)} ({rtext}× background)"
    mixed = int_list(s["mixed_lengths"])
    if st["short_state"] == "unknown":
        sv = f"{count(s['short_reads'])} reads at {qc.SHORT[0]}–{qc.SHORT[1]} nt"
    else:
        sv = f"{qc.SHORT[0]}–{qc.SHORT[1]} nt: {f2(s['short_cds_share'])} of the P-sites in the CDS"
    return [
        ("weak_frame", "Weak frame", st["weak_frame"], f"clean share {f2(s['clean_share'])}",
         f"pass ≥ {w[0]:.2f}, fail < {w[1]:.2f}", "frame"),
        ("long_reads", "Long reads", st["long_reads"], f"{pct(s['long_share'])} at ≥ {qc.LONG} nt",
         f"pass < {lr[0]:.0%}, fail ≥ {lr[1]:.0%}", "lengths"),
        ("under_digestion", "Under-digestion", st["under_digestion"],
         "too few reads at the starts" if st["under_digestion"] == "unknown" else
         f"3′ end spread {f2(s['end3_spread'])} nt",
         f"pass < {ud[0]:.1f} nt, fail ≥ {ud[1]:.1f} nt", "ends"),
        ("untemplated", "Untemplated 5′ nt", st["untemplated"], uv,
         f"fail ≥ {u[0]:.0%} and ≥ {u[1]}× background", "mismatch"),
        ("end_spike", "End spike", st["end_spike"],
         "past position 1, or at the 3′ end" if st["end_spike"] == "warn" else
         "no read sequences" if st["end_spike"] == "unknown" else "none",
         f"> {sp[0]}× background and > {sp[1]:.0%}", "mismatch"),
        ("mixed_lengths", "Mixed lengths", st["mixed_lengths"], ", ".join(f"{x} nt" for x in mixed) or "none",
         f"beside the start peak ≥ {qc.MIXED:.1f} of it", "ends"),
        ("short_state", "Short state", st["short_state"], sv,
         f"needs {qc.SHORT_READS:,} reads; warn < {qc.SHORT_CDS:.1f} in the CDS", "regions"),
        ("low_complexity", "Low complexity", st["low_complexity"],
         f"fewer than {depth} reads" if st["low_complexity"] == "unknown" else
         f"{f2(s['distinct_share_at_depth'])} distinct at {depth} reads",
         f"warn < {qc.COMPLEXITY:.2f} at {depth} reads", "depth"),
    ]


def next_steps(s, t):
    """(title, text, command or None) in order."""
    st = {f: s[f"flag_{f}"] for f in qc.FLAGS}
    lib, out = s["library"], []
    inputs = f"--gtf {s['gtf']} --fasta {s['fasta']}"
    clean = int_list(s["clean_lengths"])
    if st["untemplated"] == "fail":
        out.append(("For the next alignment: let the aligner clip the added nucleotide.",
                    f"{pct(s['mismatch5'], 0)} of the reads start with a nucleotide that does not match the transcript "
                    f"({pct(s['mismatch5_T'], 0)} of them T), and the aligner kept it in the alignment. An aligner "
                    "that soft-clips the read ends removes it: STAR with --quantTranscriptomeSAMoutput BanSingleEnd "
                    "for its transcriptome BAM (older STAR: --quantTranscriptomeBan Singleend), or bowtie2 --local. "
                    "bowtie 1 cannot soft-clip.", None))
        out.append(("For this BAM: run repair.",
                    f"repair clips the untemplated 5′ nt. Expected clean share: {f2(s['clean_share'])} → "
                    f"{f2(s['clean_share_repaired'])}.",
                    f"pitstop repair --bam {s['bam']} {inputs} --out-prefix {lib}.repaired"))
        out.append(("Inspect the repaired BAM, and take the read lengths from that report.",
                    f"Clean lengths now: {span(clean)} nt. Expected after repair: "
                    f"{span(int_list(s['clean_lengths_repaired']))} nt.",
                    f"pitstop inspect --bam {lib}.repaired.bam {inputs} --out-prefix {lib}.repaired"))
    elif s["read_lengths"]:
        lo, hi = s["read_lengths"].split("-")
        rest = [x for x in clean if not int(lo) <= x <= int(hi)]
        out.append((f"Use the clean lengths: {span(range(int(lo), int(hi) + 1))} nt.",
                    f"{f2(s['clean_share'])} of the reads are in clean lengths."
                    + (f" The other clean lengths ({span(rest)} nt) are not next to these." if rest else ""),
                    f"--read-lengths {s['read_lengths']}"))
    else:
        best = "" if isna(s["best_length"]) else (f" The best top-phase share is "
                                                   f"{f2(s['best_top_frame'])} ({int(s['best_length'])} nt).")
        out.append(("No read length has a clean frame.",
                    f"No length reaches a top-phase share of {qc.CLEAN}.{best} The steps that use the frame (orfs, score, "
                    "starts) will be weak on this library. The counts per CDS (quant) still work.", None))
    if st["untemplated"] == "info":
        out.append(("No repair needed.", f"{pct(s['mismatch5'])} of the reads have a 5′ mismatch, and the aligner "
                    f"clipped it ({pct(s['clip5_share'])} have a 5′ soft clip). pitstop removes a 5′ clip.", None))
    elif st["untemplated"] == "pass":
        out.append(("No repair needed.", f"{pct(s['mismatch5'])} of the reads have a 5′ mismatch; the background is "
                    f"{pct(s['mismatch5_background'])}.", None))
    if st["under_digestion"] in ("warn", "fail"):
        out.append(("repair does not help with under-digestion.",
                    f"The reads with the same 5′ end at the starts end at many places at the 3′ end (spread "
                    f"{f2(s['end3_spread'])} nt): the nuclease did not trim them to the ribosome. For the next "
                    "libraries, digest more. Pair this library only with libraries of the same digestion (compare "
                    "their read lengths in MultiQC).", None))
    if st["short_state"] == "warn":
        out.append((f"Leave out the {qc.SHORT[0]}–{qc.SHORT[1]} nt reads.",
                    f"Only {f2(s['short_cds_share'])} of their P-sites fall in the CDS. They are RNA fragments, not "
                    "ribosomes with an open A site.", None))
    mixed = int_list(s["mixed_lengths"])
    if mixed:
        out.append((f"Mixed lengths: {span(mixed)} nt. Check their start anchor in `pitstop starts`.",
                    "At the annotated starts, the nucleotide beside the 5′-end peak of these lengths holds "
                    f"{qc.MIXED} of the peak's reads or more: two groups of 5′ ends, and one offset puts the smaller "
                    "group 1 nt off (see the ends at the starts and stops). `pitstop starts` anchors such a length "
                    "on both nucleotides, so it can be a start length. Their frame in the CDS decides if they are "
                    "clean, as for the other lengths."
                    + (" The added 5′ nucleotide can cause this: inspect the repaired BAM."
                       if st["untemplated"] == "fail" else ""), None))
    if st["low_complexity"] == "warn":
        out.append(("Check the duplicate removal of the preprocessing.",
                    f"In a subsample of {count(qc.COMPLEXITY_DEPTH)} reads, only {f2(s['distinct_share_at_depth'])} "
                    "are distinct footprints. PCR copies and ribosomes on the same nucleotide give the same signal; "
                    "UMIs tell them apart.", None))
    if st["end_spike"] == "warn":
        out.append(("Check the trimming of the reads.",
                    "A mismatch spike past position 1 of the 5′ end, or at the 3′ end, points to left-over UMI or "
                    "adapter nucleotides. Trim them before the alignment (a fixed trim).", None))
    if st["untemplated"] == "fail" and s["untemplated_hidden"] > HIDDEN:
        out.append(("Not fixable per read.",
                    f"About {pct(s['untemplated_hidden'], 0)} of the reads carry an added nucleotide that matches the "
                    "transcript (an estimate: a third of the seen ones). No rule can find these.", None))
    return out


def inspect_report(prefix):
    """<prefix>.inspect.html from the tables of an inspect run."""
    import matplotlib
    s, t_all, a, e, c = load(prefix)
    lib = s["library"]
    t = shown_lengths(t_all)
    files = [f"{lib}.{x}" for x in ("summary.tsv", "lengths.tsv", "anchors.tsv", "ends.tsv", "complexity.tsv")]
    a5 = e[(e["length"] == "all") & (e["end"] == 5)].iloc[0]
    a3 = e[(e["length"] == "all") & (e["end"] == 3)].iloc[0]
    with matplotlib.rc_context(RC):
        charts = dict(lengths=lengths_chart(t, lib), regions=regions_chart(t, lib), frame=frame_chart(t, lib),
                      depth=depth_chart(c, lib), start=heat_chart(a, t, "start", lib),
                      stop=heat_chart(a, t, "stop", lib), mm5=mismatch_chart(a5, 5, lib),
                      mm3=mismatch_chart(a3, 3, lib))
    flags = flag_rows(s, t_all)
    st = {f: s[f"flag_{f}"] for f in qc.FLAGS}

    side = (f"<div>{LOCKUP}</div>\n<div><div class=\"kind\">inspect report</div><div class=\"lib\">{esc(lib)}</div>"
            "</div>\n<div><h3>Flags</h3><ul class=\"flaglist\">"
            + "".join(f'<li><a href="#{sec}">{icon(state)}<span>{esc(name)}</span></a></li>'
                      for _, name, state, _, _, sec in flags)
            + "</ul></div>\n<div><h3>Sections</h3><ul class=\"toc\">"
            + "".join(f'<li><a href="#{sid}">{esc(title)}</a></li>' for sid, title in SECTIONS)
            + "</ul></div>\n<div class=\"files\"><h3>Tables</h3>"
            + "".join(f"<code>{esc(f)}</code>" for f in files) + "</div>\n")

    cmd = f"pitstop inspect --bam {s['bam']} --gtf {s['gtf']} --fasta {s['fasta']} --out-prefix {lib}"
    head = ('<header><div class="eyebrow">pitstop inspect · one library</div>'
            f"<h1>{esc(lib)}</h1><dl class=\"meta\">"
            f"<dt>Command</dt><dd><code>{esc(cmd)}</code></dd>"
            f"<dt>pitstop</dt><dd>{esc(s['pitstop'])}</dd>"
            f"<dt>Reads</dt><dd>{count(s['reads'])} reads with an alignment that counts; all of them are used</dd>"
            "</dl></header>\n")

    steps = "".join(f'<li><div class="t">{esc(title)}</div><p>{esc(text)}</p>{command(cmd) if cmd else ""}</li>'
                    for title, text, cmd in next_steps(s, t_all))
    rows = pd.DataFrame(flags, columns=["Id", "Flag", "State", "Value", "Cut", "Section"])
    flag_html = ('<table class="flags"><thead><tr><th>Flag</th><th>State</th><th>Value</th><th>Cut</th><th></th>'
                 "</tr></thead><tbody>" + "".join(
                     f"<tr><td>{esc(r.Flag)}</td><td>{state_cell(r.State)}</td><td>{esc(r.Value)}</td>"
                     f'<td>{esc(r.Cut)}</td><td><a href="#{r.Section}">details</a> · '
                     f'<a href="{DOCS}{r.Id.replace("_", "-")}/">docs</a></td></tr>'
                     for r in rows.itertuples()) + "</tbody></table>")
    verdict = section("verdict", "Verdict and next steps", (
        f'<div class="vgrid"><div class="panel nextbox"><h3>Next steps</h3><ol class="steps">{steps}</ol></div>'
        f'<div class="panel scroll">{flag_html}<p class="note">The cuts were set on libraries with known problems and '
        "without. The docs page of each flag gives the reasons, and what to do. A fail does not change the exit "
        "code.</p></div></div>\n"), fold=False)

    counts = pd.DataFrame([
        ("Forward reads in the BAM", count(s["reads_forward"])),
        ("Alignments dropped for their CIGAR (indels, splices, hard clips)", count(s["alignments_dropped_cigar"])),
        ("Reads with an acceptable alignment", count(s["reads_cigar_ok"])),
        ("BAM references not in the GTF", count(s["refs_not_in_gtf"])),
        ("Alignments to them (dropped)", count(s["alignments_dropped_ref_not_in_gtf"])),
        ("Reads that count", count(s["reads"])),
        ("… on a transcript with a CDS", count(s["reads_on_cds_transcript"])),
    ], columns=["Reads", "Count"])
    nh = pd.DataFrame([(name, count(s[f"reads_nh{k}"])) for k, name in
                       (("1", "1"), ("2", "2"), ("3", "3"), ("4plus", "4 or more"))], columns=["Alignments", "Reads"])
    reads = section("reads", "The reads that pitstop reads", (
        '<p class="lead">The reads that <code>quant</code>, <code>orfs</code> and <code>starts</code> count: forward '
        "alignments whose CIGAR is matches plus soft clips at the ends, on transcripts of the GTF. A 5′ soft clip is "
        "not part of the footprint. A read with several alignments is split over them (over its alignments to "
        "transcripts with a CDS, if it has some), so each read counts once.</p>\n"
        f'<div class="figrow"><div class="panel scroll">{html_table(counts, ["Count"])}</div>'
        f'<div class="panel scroll">{html_table(nh, ["Reads"])}</div></div>\n'))

    cx = st["low_complexity"]
    d, d_at, new = s["distinct_share"], s["distinct_share_at_depth"], s["new_per_read"]
    depth_n = count(qc.COMPLEXITY_DEPTH)
    depth_find = [(cx, f"Fewer than {depth_n} reads: the distinct share at that depth is not known. Of all reads, "
                       f"{f2(d)} are distinct footprints." if cx == "unknown" else
                   f"In a subsample of {depth_n} reads, only {f2(d_at)} are distinct footprints. Check the duplicate "
                   "removal of the preprocessing." if cx == "warn" else
                   f"In a subsample of {depth_n} reads, {f2(d_at)} are distinct footprints."),
                  ("info", f"The curve still rises: {pct(new, 0)} of the last reads were new footprints. More "
                           "sequencing would add positions." if new >= RISING else
                   f"The curve flattens: {pct(new)} of the last reads were new footprints.")]
    full = c.iloc[-1]
    depth = section("depth", "Depth and complexity", (
        findings(depth_find)
        + tiles([("Reads", count(full["reads"]), "all of them"),
                 ("Distinct footprints", count(full["distinct"]), "transcript, 5′ end, length"),
                 ("Distinct share", f2(d), "all reads"),
                 (f"At {depth_n} reads", "–" if isna(d_at) else f2(d_at), f"warn < {qc.COMPLEXITY:.2f}"),
                 ("New per read", pct(new), "the last 10% of the reads")])
        + fig_block("Distinct footprints against the reads", charts["depth"],
                    "Fixed subsamples of the reads (crc32 of the read name). A curve that still rises at full depth "
                    "gains from more sequencing. A curve that is flat at a low distinct share means the library is "
                    "exhausted. Highly expressed genes give copies too, so a distinct share well below 1 is normal. "
                    "The distinct share falls with the depth, so the flag compares libraries at one depth: the "
                    f"expected distinct share of all subsamples of {depth_n} reads.",
                    legend([(None, "var(--ink)", "distinct footprints"), ("line", None, "every read distinct")]),
                    tsv_view(c, f"{lib}.complexity.tsv"))))

    sh = st["short_state"]
    short_line = (f"Too few reads at {qc.SHORT[0]}–{qc.SHORT[1]} nt ({count(s['short_reads'])}) to test the short "
                  "state." if sh == "unknown" else
                  f"{pct(s['short_share'])} of the reads are {qc.SHORT_MAX} nt or shorter, but most of the "
                  f"{qc.SHORT[0]}–{qc.SHORT[1]} nt P-sites fall outside the CDS: fragments, not the short state "
                  "(see the regions)." if sh == "warn" else
                  f"{pct(s['short_share'])} of the reads are {qc.SHORT_MAX} nt or shorter, and "
                  f"{f2(s['short_cds_share'])} of the {qc.SHORT[0]}–{qc.SHORT[1]} nt P-sites are in the CDS: ribosomes "
                  "with an open A site. Do not repair them, and do not count them as 28 nt reads.")
    disome = (f"No reads are {qc.DISOME[0]} nt or longer (the longest: {int(s['max_length'])} nt): the size "
              "selection removed disomes, so they cannot be measured here." if s["max_length"] < qc.DISOME[0] else
              f"{pct(s['disome_share'], 2)} of the reads are {qc.DISOME[0]}–{qc.DISOME[1]} nt (disomes). Expected in "
              "a disome library; in a monosome library they show how wide the size selection was.")
    shown_cls = [k for k, _ in CLASSES if (t["class"] == k).any()]
    leg = legend([(None, "var(--k-clean)", "clean"), (None, "var(--k-short)", "short state"),
                  (None, "var(--k-noframe)", "no frame"),
                  ("few", None, "too few reads")][i] for i, (k, _) in enumerate(CLASSES) if k in shown_cls)
    lt = t_all[["length", "reads", "share", "class", "mixed", "top_frame", "start_peak", "offset"]]
    lengths = section("lengths", "Read lengths", (
        findings([(st["long_reads"], f"{pct(s['long_share'])} of the reads are {qc.LONG} nt or longer. The commonest "
                                     f"length is {int(s['common_length'])} nt."),
                  (sh, short_line), ("info", disome)])
        + tiles([("Commonest length", f"{int(s['common_length'])} nt", ""),
                 (f"≥ {qc.LONG} nt", pct(s["long_share"]), f"pass < {qc.LONG_READS[0]:.0%}"),
                 (f"≤ {qc.SHORT_MAX} nt", pct(s["short_share"]), ""),
                 ("Clean share", f2(s["clean_share"]), f"lengths {span(int_list(s['clean_lengths']))}")])
        + fig_block("Reads per footprint length, by class", charts["lengths"],
                    f"Clean: top-phase share {qc.CLEAN} or more in the CDS interior. Short state: "
                    f"{qc.SHORT[0]}–{qc.SHORT[1]} nt, "
                    f"the 5′ end at nt {nt(qc.SHORT_PEAK)}, and {qc.SHORT_CDS} of the P-sites or more in the CDS. No "
                    f"frame: none of these. Too few reads: fewer than {qc.MIN_FRAME} reads in the CDS interior, or less "
                    f"than {qc.MIN_SHARE:.1%} of the reads. The "
                    f"figure shows the lengths from the first to the last with {SHOWN:.1%} of the reads; the table has all.",
                    leg, tsv_view(lt, f"{lib}.lengths.tsv"))))

    cds = s["cds_share"]
    reg_find = [("pass" if cds >= CDS_LOW else "warn",
                 f"{f2(cds)} of the P-sites fall in the CDS."
                 + ("" if cds >= CDS_LOW else f" A CDS share below {CDS_LOW} means contamination: RNA fragments or "
                                             "other RNPs."))]
    if sh == "warn":
        reg_find.append(("warn", f"{qc.SHORT[0]}–{qc.SHORT[1]} nt: only {f2(s['short_cds_share'])} of the P-sites "
                                 "fall in the CDS. Leave these lengths out."))
    elif sh == "info":
        reg_find.append(("pass", f"{qc.SHORT[0]}–{qc.SHORT[1]} nt: {f2(s['short_cds_share'])} of the P-sites fall in "
                                 "the CDS, as for ribosomes."))
    ld = s["length_distance"]
    reg_find.append(("info", "No reads on transcripts without a CDS." if isna(ld) else
                     f"The length distance between the reads on CDSs and on transcripts without a CDS is {f2(ld)}. "
                     + ("The reads without a CDS have another length profile: most are not ribosomes."
                        if ld > DISTANCE else "The two profiles are alike.")))
    rt = t_all[["length", "reads", "offset", "5UTR", "CDS", "3UTR", "transcript"]]
    regions = section("regions", "Regions per length", (
        findings(reg_find)
        + tiles([("CDS share", f2(cds), "all P-sites"),
                 (f"{qc.SHORT[0]}–{qc.SHORT[1]} nt in the CDS",
                  "–" if sh == "unknown" else f2(s["short_cds_share"]),
                  "too few reads" if sh == "unknown" else f"warn < {qc.SHORT_CDS}"),
                 ("Length distance", f2(ld), "CDS against no CDS")])
        + fig_block("Where the P-sites of each length fall", charts["regions"],
                    "The P-site of a length is its offset from the 5′ end: the distance of its commonest 5′ end from "
                    f"the start codon, or {qc.DEFAULT_OFFSET} nt with fewer than {qc.MIN_ANCHOR} reads there. The "
                    "3′ UTR starts at the stop codon.",
                    legend([(None, "var(--r-5utr)", "5′ UTR"), (None, "var(--r-cds)", "CDS"),
                            (None, "var(--r-3utr)", "3′ UTR"), (None, "var(--r-nocds)", "no CDS")]),
                    tsv_view(rt, f"{lib}.lengths.tsv"))))

    fr_find = [(st["weak_frame"], f"Clean share {f2(s['clean_share'])}: the clean lengths are "
                                  f"{span(int_list(s['clean_lengths']))} nt."
                + ("" if isna(s["best_length"]) else
                   f" The best top-phase share is {f2(s['best_top_frame'])} at {int(s['best_length'])} nt."))]
    if abs(s["clean_share_repaired"] - s["clean_share"]) >= CHANGED:
        fr_find.append(("info", f"After repair the clean share would be {f2(s['clean_share_repaired'])} "
                                f"({span(int_list(s['clean_lengths_repaired']))} nt)."))
    ft = t_all[["length", "interior", "phase0", "phase1", "phase2", "top_phase", "top_frame", "class",
                "top_frame_repaired", "class_repaired"]]
    frame = section("frame", "Frame per length", (
        findings(fr_find)
        + fig_block("5′ phase of the reads in the CDS interior", charts["frame"],
                    "Each bar is one length; its bottom segment is the top phase, so a bottom segment that crosses "
                    f"the dashed line ({qc.CLEAN}) is a clean length. Faded: too few reads. A ring marks the top-phase "
                    "share after repair, where it differs. The CDS interior: the read is at least "
                    f"{qc.INTERIOR} nt inside the CDS at both ends.",
                    legend([(None, "var(--f0)", "phase 0"), (None, "var(--f1)", "phase 1"),
                            (None, "var(--f2)", "phase 2"), ("ring", None, "top-phase share after repair"),
                            ("line", None, f"clean: {qc.CLEAN}")]),
                    tsv_view(ft, f"{lib}.lengths.tsv"))))

    ends = section("ends", "Ends at the annotated starts and stops", (
        findings(end_lines(s, t_all, st))
        + '<div class="figrow">'
        + fig_block("5′ ends around the start codon", charts["start"],
                    "The read's 5′ end; nt 0 is the first nt of the start codon. Each row is one length, and its "
                    "cells add up to 1.", small=True)
        + fig_block("3′ ends around the stop codon", charts["stop"],
                    "The read's 3′ end; nt 0 is the first nt of the stop codon. Extra length at the 3′ end moves a "
                    "row to the right.", small=True)
        + "</div>\n"
        + legend([(None, f"var(--h{i})", f"{'< ' + str(HEAT_BINS[0]) if i == 0 else '≥ ' + str(HEAT_BINS[i - 1])}")
                  for i in range(8)])
        + tsv_view(t_all[["length", "start_reads", "start_peak", "start_peak_share", "start_beside", "stop_reads",
                          "stop_peak", "stop_peak_share", "stop_beside", "class", "mixed"]], f"{lib}.lengths.tsv")
        + tsv_view(a[a["length"].isin(t["length"])], f"{lib}.anchors.tsv")))

    mismatch = section("mismatch", "Read ends against the transcript", mismatch_section(s, e, a5, charts, st, lib))

    scope = section("scope", "Not in this report", (
        '<ul class="scope">'
        "<li><strong>Raw reads:</strong> base quality, adapter content, over-represented sequences. Use FastQC on the "
        "FASTQ and the log of the adapter trimmer.</li>"
        "<li><strong>rRNA and unaligned reads:</strong> the BAM holds only the reads that aligned to transcripts. Use "
        "the log of the rRNA filter and of the aligner.</li>"
        "<li><strong>Codon dwell times:</strong> they need P-site offsets, and they describe biology, not the "
        "library.</li>"
        "<li><strong>FLOSS per transcript:</strong> a score for single ORFs, not for a library. This report gives the "
        "library view (regions per length).</li>"
        "<li><strong>Several libraries:</strong> MultiQC reads the <code>_mqc</code> files of each run and shows the "
        "libraries side by side. <code>starts</code> checks the pairs of initiation and elongation libraries.</li>"
        "</ul>\n"), fold=False)
    foot = ("<footer><div>Tables: " + ", ".join(f"<code>{esc(f)}</code>" for f in files) + "</div><div>Each number "
            "in this report is also in a table. A rerun with the same inputs gives a byte-identical report.</div>"
            "</footer>\n")
    out = page(f"{lib} · pitstop inspect", side, head + verdict + reads + depth + lengths + regions + frame + ends
               + mismatch + scope + foot)
    path = f"{prefix}.inspect.html"
    with open(path, "w", encoding="utf-8") as f:
        f.write(out)
    log.info("wrote %s", path)
    multiqc(prefix, s, t_all)
    return path


def end_lines(s, t, st):
    """The findings of the section on the ends at the starts and stops."""
    out = []
    idx = t.set_index("length")
    core = idx.loc[int(s["common_length"])]
    long = idx[(idx.index >= qc.LONG) & (idx["start_reads"] > 0) & idx["start_peak"].notna()]
    long = long[long["start_peak_share"] * long["start_reads"] >= qc.MIN_ANCHOR]
    ud = st["under_digestion"]
    if ud != "unknown":
        out.append((ud, f"At the starts, the reads whose 5′ end is at nt {nt(s['end3_peak'])} (longer than "
                        f"{qc.SHORT_MAX} nt) have a 3′ end spread of {f2(s['end3_spread'])} nt: the standard deviation "
                        "of their lengths. "
                        + ("Their 3′ ends are at many places: under-digestion. No rule per read fixes this."
                           if ud != "pass" else "Their 3′ ends are trimmed to a few places.")))
    if st["long_reads"] != "pass" and len(long) and not isna(core["start_peak"]):
        moved = (long["start_peak"] - core["start_peak"]).abs()
        peaks = ", ".join(nt(x) for x in sorted(set(long["start_peak"]), reverse=True))
        three = bool((moved <= 1).all())
        out.append((st["long_reads"],
                    f"At the starts, the 5′ ends of the reads of {qc.LONG} nt or more peak at nt {peaks}; those of the "
                    f"commonest length ({int(s['common_length'])} nt) at nt {nt(core['start_peak'])}. "
                    + ("So the extra length is at the 3′ end." if three else
                       "So part of the extra length is at the 5′ end.")))
    if st["untemplated"] == "fail":
        out.append(("warn", "Part of each length has its 5′ end 1 nt upstream: the added 5′ nucleotide, not a second "
                            "group of footprints. repair clips the untemplated 5′ nt and moves these reads back."))
    clean = idx[idx["class"] == "clean"]
    if len(clean):
        peaks = sorted({nt(x) for x in clean["start_peak"].dropna()})
        out.append(("pass" if len(peaks) == 1 else "info",
                    f"At the starts, the 5′ ends of the clean lengths peak at nt {', '.join(peaks)}."))
    mixed = idx[idx["mixed"]]
    if len(mixed):
        each = "; ".join(f"{x} nt: {f2(r['start_beside'])} beside nt {nt(r['start_peak'])}" for x, r in mixed.iterrows())
        out.append(("info", f"Mixed lengths ({each}): the nucleotide beside the 5′-end peak holds this share of the "
                            "peak's reads. One offset per length uses the larger group."))
    return out or [("info", "Too few reads at the annotated starts to read the ends.")]


def mismatch_section(s, e, a5, charts, st, lib):
    u = st["untemplated"]
    m5, bg5 = s["mismatch5"], s["mismatch5_background"]
    if u == "fail":
        line = (f"The 5′ mismatch is {pct(m5)} against a background of {pct(bg5)}; {pct(s['mismatch5_T'], 0)} of these "
                f"nucleotides are T. The aligner did not clip them (5′ clip {pct(s['clip5_share'])}). An aligner "
                "that soft-clips the read ends removes them in the next alignment; for this BAM, repair clips them "
                "(Next steps).")
    elif u == "info":
        line = (f"The 5′ mismatch is {pct(m5)}, and the aligner clipped it (5′ clip {pct(s['clip5_share'])}). pitstop "
                "removes a 5′ clip: nothing to do.")
    elif u == "unknown":
        line = "The BAM has no read sequences, so the read ends cannot be compared with the transcript."
    else:
        line = f"The 5′ mismatch ({pct(m5)}) is near the background ({pct(bg5)}). No added nucleotide to repair."
    sp = st["end_spike"]
    items = [(u, line), (sp, "A mismatch spike past position 1 of the 5′ end, or at the 3′ end: check the trimming "
                             "(a fixed trim before the alignment)." if sp == "warn" else
                         "No mismatch spike past position 1, and none at the 3′ end: no fixed trim needed."
                         if sp == "pass" else "No read sequences.")]
    b = s["bias_max"]
    if not isna(b):
        strong = abs(np.log2(b)) >= np.log2(BIAS)
        items.append(("info", f"The strongest end preference is {b:.2f}× ({esc(s['bias_max_at'])} at the cut)"
                      + (". The counts carry this bias. For the next libraries, check the nuclease and the ligation."
                         if strong else ": weak, no action.")))
    bases = "".join(
        f'<span class="b">{x}</span><span class="track"><span class="fill" style="width:'
        f'{0 if isna(a5[f"mm0_{x}"]) else 100 * a5[f"mm0_{x}"]:.1f}%;display:block"></span></span>'
        f'<span class="r">{pct(a5[f"mm0_{x}"], 0)}</span>' for x in qc.BASES)
    all_rows = e[e["length"] == "all"]
    bias = pd.DataFrame([[f"{'5′' if r['end'] == 5 else '3′'} cut"] + [f2(r[f"bias_{x}"]) for x in qc.BASES]
                         for _, r in all_rows.iterrows()], columns=["End"] + list(qc.BASES))
    ecols = ["length", "end", "reads", "clip1", "clip2", "clip3plus"] + [f"mm{i}" for i in range(qc.K)] \
        + [f"mm0_{x}" for x in qc.BASES] + [f"cut_{x}" for x in qc.BASES] + [f"bias_{x}" for x in qc.BASES]
    return (findings(items)
            + tiles([("5′ mismatch", pct(m5), "position 0, the clip included"),
                     ("Background", pct(bg5), f"positions {qc.BACKGROUND[0]}–{qc.BACKGROUND[1] - 1}"),
                     ("5′ soft clip", pct(s["clip5_share"]), "clipped by the aligner"),
                     ("Hidden added nt", pct(s["untemplated_hidden"]), "an estimate")])
            + '<div class="figrow">'
            + fig_block("Mismatch rate from the 5′ end", charts["mm5"],
                        "Position 0 is the first nucleotide of the read, the soft clip included. Dashed: the "
                        "background.", small=True)
            + fig_block("Mismatch rate from the 3′ end", charts["mm3"],
                        "A spike here, or over more than 2 positions at the 5′ end, points to left-over UMI or adapter "
                        "nucleotides.", small=True)
            + "</div>\n"
            + f'<figure><div class="ft">The read\'s nucleotide at a mismatched 5′ position 0</div>'
              f'<div class="bases">{bases}</div><figcaption>Reverse transcriptase adds mostly T at the 5′ end of the '
              "read. An even split points to something else.</figcaption></figure>\n"
            + '<figure><div class="ft">The transcript nucleotide at each cut, against one codon inside</div>'
              f'<div class="scroll">{html_table(bias, list(qc.BASES), "bias")}</div><figcaption>The first and the '
              "last nucleotide of the footprint (after the 5′ clip), as a share, over the share of the nucleotide 3 nt "
              "inside the footprint (the same codon position, not at a cut). A ratio far from 1 is a preference of "
              "the nuclease or the ligation. It changes the counts; the advice is for the next libraries."
              f"</figcaption>{tsv_view(e[ecols], f'{lib}.ends.tsv')}</figure>\n")


def multiqc(prefix, s, t):
    """MultiQC custom content (found by its _mqc names, no MultiQC config): a row of the general statistics, and two
    line graphs (the read lengths, the top-phase share per length). The files of all libraries go into one section each."""
    lib = s["library"]
    st = [s[f"flag_{f}"] for f in qc.FLAGS]
    r = lambda v: None if isna(v) else round(float(v), 4)
    general = {
        "id": "pitstop_inspect_general", "plot_type": "generalstats",
        "pconfig": [
            {"clean_share": {"title": "Clean share", "description": "pitstop inspect: reads in clean lengths",
                             "min": 0, "max": 1, "format": "{:,.2f}", "scale": "RdYlGn"}},
            {"clean_share_repaired": {"title": "Clean (repair)", "description": "pitstop inspect: the clean share "
                                      "after repair", "min": 0, "max": 1, "format": "{:,.2f}",
                                      "scale": "RdYlGn"}},
            {"mismatch5": {"title": "5' mismatch", "description": "pitstop inspect: the mismatch rate at the first nt "
                           "of the read", "min": 0, "max": 100, "suffix": "%", "format": "{:,.1f}", "scale": "OrRd"}},
            {"long_share": {"title": f">= {qc.LONG} nt", "description": f"pitstop inspect: reads of {qc.LONG} nt or "
                            "more", "min": 0, "max": 100, "suffix": "%", "format": "{:,.1f}", "scale": "OrRd"}},
            {"distinct_share_at_depth": {"title": f"Distinct ({qc.COMPLEXITY_DEPTH / 1e6:g}M)",
                                         "description": "pitstop inspect: distinct footprints per read in a subsample "
                                         f"of {qc.COMPLEXITY_DEPTH:,} reads", "min": 0, "max": 1,
                                         "format": "{:,.2f}", "scale": "RdYlGn"}},
            {"flags_fail": {"title": "Fails", "description": "pitstop inspect: flags that fail", "min": 0,
                            "format": "{:,.0f}", "scale": "Reds"}},
            {"flags_warn": {"title": "Warns", "description": "pitstop inspect: flags that warn", "min": 0,
                            "format": "{:,.0f}", "scale": "Oranges"}},
        ],
        "data": {lib: {"clean_share": r(s["clean_share"]), "clean_share_repaired": r(s["clean_share_repaired"]),
                       "mismatch5": r(100 * s["mismatch5"]) if not isna(s["mismatch5"]) else None,
                       "long_share": r(100 * s["long_share"]),
                       "distinct_share_at_depth": r(s["distinct_share_at_depth"]),
                       "flags_fail": st.count("fail"), "flags_warn": st.count("warn")}},
    }
    lengths = {
        "id": "pitstop_inspect_lengths", "section_name": "pitstop inspect: read lengths",
        "description": "The share of the reads per footprint length (pitstop inspect, lengths.tsv).",
        "plot_type": "linegraph",
        "pconfig": {"id": "pitstop_inspect_lengths_plot", "title": "pitstop inspect: read lengths",
                    "xlab": "footprint length (nt)", "ylab": "share of the reads", "ymin": 0},
        "data": {lib: [[int(x), r(v)] for x, v in zip(t["length"], t["share"])]},
    }
    framed = t[t["class"] != "few reads"]
    frame = {
        "id": "pitstop_inspect_top_frame", "section_name": "pitstop inspect: top-phase share per length",
        "description": "The share of the CDS-interior reads whose 5' end is in the commonest phase, per footprint "
                       "length (1/3: no frame; 1: all reads in one phase), for the lengths with "
                       f"{qc.MIN_FRAME} reads or more in the CDS interior (pitstop inspect, lengths.tsv: top_frame). "
                       f"Clean: {qc.CLEAN} or more.",
        "plot_type": "linegraph",
        "pconfig": {"id": "pitstop_inspect_top_frame_plot", "title": "pitstop inspect: top-phase share per length",
                    "xlab": "footprint length (nt)", "ylab": "top-phase share", "ymin": 0, "ymax": 1},
        "data": {lib: [[int(x), r(v)] for x, v in zip(framed["length"], framed["top_frame"])]},
    }
    for name, obj in (("general", general), ("lengths", lengths), ("frame", frame)):
        with open(f"{prefix}.inspect_{name}_mqc.json", "w") as f:
            json.dump(obj, f, indent=1)
            f.write("\n")

"""The report of `starts` (<prefix>.starts.html), from the tables of its run.

A run pools pairs (a row of the library sheet: an initiation library and its elongation library) by pool, so the
start lengths and the calls belong to the run, and the QC of the inputs belongs to each pair: the same digestion (the
read-length distance) and the drug effect (the start peak). The page: the verdict of the run, the pairs (a
table, a figure, and a tab per pair with the two libraries side by side), the start lengths, harringtonine at the
annotated starts per pool, the calls.
"""
import logging
import os

import numpy as np
import pandas as pd

from .. import starts as sx
from ..logo import LOCKUP
from .base import (C, HEAT_BINS, RC, STATE, count, esc, f2, fig_block, findings, html_table, icon, int_list, isna,
                   legend, nt, page, pct, section, span, state_cell, svg, tiles, tsv_view)

log = logging.getLogger(__name__)

S_SECTIONS = (("verdict", "Verdict and next steps"), ("calls", "Calls"), ("pairs", "Pairs"), ("lengths", "Start lengths"),
              ("pair", "Pair detail"), ("scope", "Not in this report"))
REGION_VAR = (("5UTR", "r-5utr", "5′ UTR"), ("CDS", "r-cds", "CDS"), ("3UTR", "r-3utr", "3′ UTR"),
              ("transcript", "r-nocds", "no CDS"))
CALL_CLASSES = (("annotated", "Annotated starts"), ("5UTR", "5′ UTR"), ("CDS", "CDS"), ("3UTR", "3′ UTR"),
                ("transcript", "No CDS"))
ORF_TYPES = ("uORF", "uoORF", "CDS", "extension", "truncation", "internal", "dORF", "other")
VISIBLE = 0.004      # the profile bars draw from this share of the row's highest bar


# ---------- the tables ----------

def load_starts(prefix):
    """The tables of a starts run: the stats (a dict of strings), lengths, profiles, pairs, and libraries, kernel and
    start codons (None when no read length passed the QC: the run stopped before them)."""
    stats = (pd.read_csv(f"{prefix}.starts_stats.tsv", sep="\t", dtype=str, keep_default_na=False)
             .set_index("stat")["value"].to_dict())
    t = pd.read_csv(f"{prefix}.lengths.tsv", sep="\t", dtype={"library": str, "pool": str, "qc": str})
    prof = pd.read_csv(f"{prefix}.profiles.tsv", sep="\t", dtype={"library": str})
    pairs = pd.read_csv(f"{prefix}.pairs.tsv", sep="\t", dtype={"pool": str, "harringtonine": str, "elongation": str})
    done = bool(stats.get("start_lengths"))
    read = lambda name, **kw: pd.read_csv(f"{prefix}.{name}.tsv", sep="\t", **kw) if done else None
    return dict(stats=stats, t=t, prof=prof, pairs=pairs, lt=read("libraries", dtype={"library": str, "pool": str}),
                kern=read("kernel", dtype={"pool": str}), codons=read("start_codons", keep_default_na=False))


def library_order(t):
    """One row per library: by pool, harringtonine before elongation, then the order of the table; its file name."""
    libs = t.drop_duplicates("library")[["library", "pool", "kind"]].reset_index(drop=True)
    libs["elong"] = libs["kind"] != "harringtonine"
    libs = libs.sort_values(["pool", "elong"], kind="stable").reset_index(drop=True)
    libs["name"] = libs["library"].map(os.path.basename)
    return libs


def pair_order(pairs, lt):
    """The pairs by pool (then the order of the sheet): a label, the file names, and the states of the two checks of
    a pair: the read-length distance, and the drug effect (unknown without libraries.tsv)."""
    p = pairs.sort_values("pool", kind="stable").reset_index(drop=True)
    p["n"] = p.groupby("pool").cumcount() + 1
    p["label"] = [f"{pl} · pair {n}" for pl, n in zip(p["pool"], p["n"])]
    p["h_name"], p["e_name"] = p["harringtonine"].map(os.path.basename), p["elongation"].map(os.path.basename)
    p["digestion"] = np.where(p["distance"] >= sx.PAIR_DISTANCE, "warn", "pass")
    drug = []
    for h, e in zip(p["harringtonine"], p["elongation"]):
        if lt is None:
            drug.append("unknown")
            continue
        a, b = lt.set_index("library").loc[h], lt.set_index("library").loc[e]
        if isna(a["start_enrichment"]) or isna(b["start_enrichment"]):
            drug.append("unknown")
        else:
            drug.append("warn" if a["start_enrichment"] <= b["start_enrichment"] else "pass")
    p["drug"] = drug
    return p


def names(df, *columns):
    """df with the library prefixes of `columns` cut to their file names: a report has no paths."""
    return df.assign(**{c: df[c].map(os.path.basename) for c in columns})


def length_window(stats):
    lo, hi = stats["length_window"].split("-")
    return int(lo), int(hi)


def window_lengths(t, stats):
    """The read lengths the report draws: those in the length window, and the start lengths."""
    lo, hi = length_window(stats)
    keep = t["length"].between(lo, hi) | t["start_length"]
    return sorted(int(x) for x in t.loc[keep, "length"].unique()) or sorted(int(x) for x in t["length"].unique())


def moved_lengths(t, chosen, library=None, pair=False):
    """The start lengths that the calibration moves in a harringtonine library (or in `library`): [(length, shift)].
    pair: the pair lengths only, else the others."""
    m = t[(t["kind"] == "harringtonine") & t["length"].isin(chosen) & (t["shift"] != 0) & (t["pair"] == pair)]
    if library is not None:
        m = m[m["library"] == library]
    return sorted({(int(x), int(v)) for x, v in zip(m["length"], m["shift"])})


# ---------- numbers as text ----------

def compact(v):
    """A tick label: 950, 12k, 1.2M."""
    return (f"{v / 1e6:.1f}M" if v >= 1e6 else f"{v / 1e3:.0f}k" if v >= 1e4 else f"{v / 1e3:.1f}k" if v >= 1e3
            else f"{v:.0f}")


def dist(x):
    """A read-length distance: 2 decimals, or 4 where 2 would round onto the cut."""
    v = f"{x:.2f}"
    return f"{x:.4f}" if float(v) == sx.PAIR_DISTANCE and x != sx.PAIR_DISTANCE else v


def fold(x):
    """An enrichment: 34.8×."""
    return "–" if isna(x) else f"{x:.1f}×"


# ---------- figures ----------

def bars(ax, x, y, colours, width=0.85, alpha=1.0, base=None):
    """Bars from base (default 0) to base + y (also below), each a stroked line of the bar's width, one path per
    colour: a page with thousands of bars stays small. Set the x limits first: the line width comes from them."""
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    x0, x1 = ax.get_xlim()
    pt = ax.get_position().width * ax.figure.get_figwidth() * 72 * width / (x1 - x0)
    floor = max(VISIBLE * np.abs(y).max(), 1e-300)   # lower bars are not visible
    base = np.zeros(len(x)) if base is None else base
    for c in dict.fromkeys(colours):
        k = [i for i in range(len(x)) if colours[i] == c and abs(y[i]) >= floor]
        if not k:
            continue
        verts = [v for i in k for v in ((x[i], base[i]), (x[i], base[i] + y[i]))]
        ax.add_patch(PathPatch(Path(verts, [Path.MOVETO, Path.LINETO] * len(k)), facecolor="none", edgecolor=c,
                               linewidth=pt, capstyle="butt", alpha=alpha))


def cells(ax, xs, rows, values, gap=0.06):
    """A heat map: a cell per (row, x), in the colour of its step of HEAT_BINS (h0-h7), one path per step with a gap
    between the cells (the background shows through): a few paths, not one per cell."""
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    step = np.searchsorted(HEAT_BINS, values, side="right")
    for b in range(len(HEAT_BINS) + 1):
        r, c = np.nonzero(step == b)
        if not len(r):
            continue
        verts, codes = [], []
        for y, x in zip(rows[r], xs[c]):
            lo, hi, top, bot = x - 0.5 + gap, x + 0.5 - gap, y - 0.5 + gap, y + 0.5 - gap
            verts += [(lo, top), (hi, top), (hi, bot), (lo, bot), (lo, top)]
            codes += [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=C[f"h{b}"], edgecolor="none", linewidth=0))


def share_bar(ax, r, parts, label_part, what, name, gid, tips):
    """One horizontal bar of shares (a row of libraries.tsv), the share of part label_part written in its segment."""
    left = 0.0
    for j, (key, var, label) in enumerate(parts):
        v = float(r[key])
        if v <= 0:
            continue
        b = ax.barh([0], [v], left=left, height=0.62, color=C[var], edgecolor=C["surface"], linewidth=0.8)
        b[0].set_gid(f"{gid}-{j}")
        tips[f"{gid}-{j}"] = f"{name}, {label}: {pct(v)} of the start-length P-sites{what}"
        if j == label_part and v >= 0.12:
            ax.text(left + v / 2, 0, pct(v, 0), ha="center", va="center", color=C["surface"], fontsize=6.5,
                    fontweight="bold")
        left += v
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def pairs_chart(p, lt, key):
    """A row per pair: the read-length distance against its cut; with libraries.tsv, the start peaks of the two
    libraries (log scale)."""
    from matplotlib.figure import Figure
    n, side = len(p), lt is not None
    h = 0.75 + 0.24 * n
    fig = Figure(figsize=(8.0, h))
    gs = fig.add_gridspec(1, 2 if side else 1, left=0.18, right=0.98,
                          top=1 - 0.3 / h, bottom=0.45 / h, wspace=0.12)
    rows = np.arange(n)
    tips = {}
    ax = fig.add_subplot(gs[0, 0])
    for i, r in p.iterrows():
        b = ax.barh([i], [r["distance"]], height=0.6, linewidth=0,
                    color=C["warn"] if r["digestion"] == "warn" else C["k-short"])
        b[0].set_gid(f"{key}-d{i}")
        tips[f"{key}-d{i}"] = f"{r['label']}: length distance {dist(r['distance'])}"
    ax.axvline(sx.PAIR_DISTANCE, color=C["ink2"], linewidth=1, linestyle=(0, (4, 3)))
    ax.set_xlim(0, max(1.0, float(p["distance"].max()) * 1.05))
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_yticks(rows, p["label"])
    ax.set_xlabel("read-length distance")
    ax.set_title("same digestion", loc="left", fontsize=7.5, color=C["ink2"], pad=4)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    if side:
        s = lt.set_index("library")
        for col, what, title in ((1, "start_enrichment", "start peak"),):
            ax = fig.add_subplot(gs[0, col])
            for i, r in p.iterrows():
                a, b = s.at[r["harringtonine"], what], s.at[r["elongation"], what]
                if isna(a) or isna(b):
                    continue
                ax.plot([b, a], [i, i], color=C["axis"], linewidth=1.5, zorder=1)
                ax.scatter([b], [i], s=30, color=C["k-short"], edgecolor=C["surface"], linewidth=0.8, zorder=2)
                pt = ax.scatter([a], [i], s=30, color=C["k-clean"], edgecolor=C["surface"], linewidth=0.8, zorder=3)
                pt.set_gid(f"{key}-{col}-{i}")
                tips[f"{key}-{col}-{i}"] = f"{r['label']}: {title} {fold(a)} initiation, {fold(b)} elongation"
            from matplotlib.ticker import FixedLocator, NullFormatter
            vals = s.loc[list(p["harringtonine"]) + list(p["elongation"]), what].dropna()
            lo_x, hi_x = (min(float(vals.min()), 1) / 1.6, max(float(vals.max()), 1) * 1.6) if len(vals) else (0.5, 2)
            ax.set_xscale("log")
            ax.set_xlim(lo_x, hi_x)
            ax.xaxis.set_major_locator(FixedLocator([v for v in (0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000)
                                                     if lo_x <= v <= hi_x]))
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.axvline(1, color=C["ink2"], linewidth=0.8, linestyle=(0, (2, 2)))
            ax.set_ylim(n - 0.5, -0.5)
            ax.set_yticks(rows, [""] * n)
            ax.set_xlabel(f"{title} (× the mean in-frame codon)")
            ax.set_title(title, loc="left", fontsize=7.5, color=C["ink2"], pad=4)
            ax.xaxis.set_major_formatter(lambda v, _: f"{v:g}")
            ax.grid(axis="x")
            ax.set_axisbelow(True)
    return svg(fig, key, "Read-length distance and start peaks per pair", tips)


PAIR_LIBS = (("harringtonine", "initiation", -0.2, 1.0), ("elongation", "elongation", 0.2, 0.5))


def lengths_pair_chart(t, r, key):
    """The read lengths of the two libraries of a pair: their shares of CDS-interior P-sites, side by side; below,
    the frames of those P-sites per length (initiation left, elongation right and faded)."""
    from matplotlib.figure import Figure
    fig = Figure(figsize=(8.0, 3.9))
    gs = fig.add_gridspec(2, 1, left=0.08, right=0.99, top=0.95, bottom=0.11, hspace=0.25)
    ax, fx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    tips = {}
    top = 1e-9
    L = sorted(t["length"].unique())
    for a in (ax, fx):
        a.set_xlim(L[0] - 0.7, L[-1] + 0.7)
        a.set_xticks(L if len(L) <= 26 else L[::2])
        a.grid(axis="y")
        a.set_axisbelow(True)
    for j, (col, kind, dx, alpha) in enumerate(PAIR_LIBS):
        q = t[t["library"] == r[col]]
        colour = C["k-clean"] if j == 0 else C["k-short"]
        for x, v in zip(q["length"], q["interior_share"]):
            b = ax.bar([x + dx], [v], width=0.38, color=colour, linewidth=0)
            b[0].set_gid(f"{key}-{j}-{x}")
            tips[f"{key}-{j}-{x}"] = f"{x} nt, {kind}: {pct(v)} of the CDS-interior P-sites"
        top = max(top, float(q["interior_share"].max()))
        fr = q[["frame0", "frame1", "frame2"]].to_numpy(dtype=float)
        fr = fr / np.maximum(fr.sum(axis=1, keepdims=True), 1)
        xs, low = q["length"].to_numpy() + dx, np.zeros(len(q))
        for f in range(3):
            bars(fx, xs, fr[:, f], [C[f"f{f}"]] * len(q), width=0.38, alpha=alpha, base=low)
            low = low + fr[:, f]
    ax.tick_params(labelbottom=False)
    ax.set_ylim(0, top * 1.12)
    ax.set_ylabel("share of P-sites")
    ax.yaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%")
    ax.set_title("read lengths", loc="left", fontsize=7.5, color=C["ink2"], pad=3)
    fx.set_ylim(0, 1)
    fx.set_yticks([0, 0.5, 1])
    fx.yaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%")
    fx.set_ylabel("share by frame")
    fx.set_xlabel("read length (nt)")
    fx.set_title("frames", loc="left", fontsize=7.5, color=C["ink2"], pad=3)
    return svg(fig, key, f"Read lengths and frames of {r['h_name']} and {r['e_name']}", tips)


def pair_profiles(t, prof, r, L):
    """The tables of a pair's profile figures: lengths.tsv per library by length, and profiles.tsv per anchor."""
    q, tab = {}, {}
    for col, _, _, _ in PAIR_LIBS:
        q[col] = t[t["library"] == r[col]].set_index("length").reindex(L)
        p = prof[prof["library"] == r[col]]
        tab[col] = {a: g.pivot_table(index="length", columns="x", values="psites", aggfunc="sum")
                    .reindex(L, fill_value=0) for a, g in p.groupby("anchor")}
    return q, tab


def panel_title(ax, length, start):
    ax.set_title(f"{length} nt" + (" · start length" if start else ""), loc="left", fontsize=8.5, pad=4,
                 fontweight="bold" if start else "normal", color=C["ink"] if start else C["muted"])


def profile_rows(t, prof, r, L, anchor, key):
    """A row per read length: the P-sites of the pair per nt at the annotated starts or at the stop codons, coloured
    by the frame of the CDS. Initiation above, elongation below (faded, hanging down): two axes with a gap between
    them, each with its own scale, its label and a tick at its highest bar."""
    from matplotlib.figure import Figure
    q, tab = pair_profiles(t, prof, r, L)
    h_, e_ = (tab[col][anchor] for col in ("harringtonine", "elongation"))
    x = h_.columns.to_numpy()
    e_ = e_.reindex(columns=x, fill_value=0)
    colours = [C[f"f{v % 3}"] for v in x]
    ticks = [v for v in x if v % 15 == 0]
    title, half, gap, row = 0.24, 0.5, 0.06, 1.64    # inches: a row is the title, two halves and the x ticks
    h = 0.08 + row * len(L) + 0.2
    fig = Figure(figsize=(8.0, h))
    left, width = 0.95 / 8.0, 6.95 / 8.0
    for i, length in enumerate(L):
        a, b = q["harringtonine"].loc[length], q["elongation"].loc[length]
        top = h - 0.08 - row * i - title
        up = fig.add_axes((left, (top - half) / h, width, half / h))
        down = fig.add_axes((left, (top - 2 * half - gap) / h, width, half / h), sharex=up)
        for ax, y, alpha, lib, kind, ink in ((up, h_, 1.0, a, "initiation", C["ink"]),
                                              (down, e_, 0.5, b, "elongation", C["ink2"])):
            y = y.loc[length].to_numpy(dtype=float)
            m = max(y.max(), 1)
            ax.set_xlim(x[0] - 1, x[-1] + 1)
            ax.axvline(0, color=C["axis"], linewidth=0.6, linestyle=(0, (2, 2)), zorder=0)
            bars(ax, x, y, colours, alpha=alpha)
            ax.set_ylim((0, m * 1.18) if ax is up else (m * 1.18, 0))
            ax.set_yticks([m], [compact(m)])
            ax.grid(axis="y", linestyle=(0, (1, 2)))
            ax.set_axisbelow(True)
            ax.tick_params(axis="y", labelcolor=ink, labelsize=7)
            ax.set_ylabel(kind, rotation=0, ha="right", va="center", color=ink, fontsize=7.5, labelpad=6)
            ax.spines["bottom" if ax is up else "top"].set_visible(True)
            ax.spines["bottom" if ax is up else "top"].set_color(C["axis"])
            if anchor == "start" and lib["at_peak"]:
                ax.text(0.995, 0.96 if ax is up else 0.04,
                        f"peak {nt(lib['start_peak'])}" + (", mixed" if lib["mixed"] else ""), ha="right",
                        va="top" if ax is up else "bottom", fontsize=7, color=C["ink2"], transform=ax.transAxes)
        up.tick_params(axis="x", bottom=False, labelbottom=False)
        down.set_xticks(ticks, [nt(v) for v in ticks])
        panel_title(up, length, bool(a["start_length"]))
        if i == len(L) - 1:
            down.set_xlabel("nt from the start codon" if anchor == "start" else "nt from the stop codon")
    what = "annotated starts" if anchor == "start" else "stop codons"
    return svg(fig, key, f"P-sites of {r['h_name']} and {r['e_name']} at the {what}, per read length")


def cds_grid(t, prof, r, L, key):
    """A panel per read length, four per row: the P-sites of the pair along the CDS (a line per library, each as a
    share of its P-sites in the CDS)."""
    from matplotlib.figure import Figure
    q, tab = pair_profiles(t, prof, r, L)
    n, cols = len(L), 4
    rows = -(-n // cols)
    h = 0.25 + 1.15 * rows
    fig = Figure(figsize=(8.0, h))
    gs = fig.add_gridspec(rows, cols, left=0.02, right=0.99, top=1 - 0.28 / h, bottom=0.3 / h, hspace=0.75,
                          wspace=0.1)
    for i, length in enumerate(L):
        ax = fig.add_subplot(gs[i // cols, i % cols])
        for col, colour, style in (("harringtonine", C["ink"], "-"), ("elongation", C["muted"], (0, (3, 2)))):
            c = tab[col].get("cds")
            y = c.loc[length].to_numpy(dtype=float) if c is not None else np.zeros(sx.CDS_BINS)
            ax.plot(np.arange(len(y)), y / max(y.sum(), 1), color=colour, linewidth=1.3, linestyle=style)
        ax.set_yticks([])
        ax.set_xticks([0, sx.CDS_BINS - 1], ["start", "stop"])
        ax.set_ylim(bottom=0)
        panel_title(ax, length, bool(q["harringtonine"].loc[length, "start_length"]))
    return svg(fig, key, f"P-sites of {r['h_name']} and {r['e_name']} along the CDS, per read length")


def length_table(t, r, L):
    """The pair per read length, as text: the share of the library's P-sites, the offset (and the calibration's
    move), the start peak (and the nt beside it, if mixed), the top-phase share inside the CDS. (The QC of a length
    is the run's, in the section Start lengths.)"""
    q = {col: t[t["library"] == r[col]].set_index("length").reindex(L) for col in ("harringtonine", "elongation")}
    cells = lambda x: [pct(x["share"]),
                       "–" if isna(x["offset"]) else f"{int(x['offset'])}" + (f" ({nt(x['shift'])})" if x["shift"]
                                                                               else ""),
                       ("–" if not x["at_peak"] else nt(x["start_peak"])
                        + (f", mixed: {nt(x['beside_nt'])} has {f2(x['beside_share'])}" if x["mixed"] else "")),
                       f2(x["top_frame"])]
    head = ("<thead><tr><th></th><th class=grp colspan=4>Initiation</th><th class=grp colspan=4>Elongation</th>"
            "</tr><tr><th class=r>Length</th>"
            + "<th class=r>P-sites</th><th class=r>Offset</th><th>Start peak</th><th class=r>Top phase</th>" * 2
            + "</tr></thead>")
    body = []
    for length in L:
        a, b = q["harringtonine"].loc[length], q["elongation"].loc[length]
        start = bool(a["start_length"])
        tds = [f"{length} nt"] + cells(a) + cells(b)
        body.append(f'<tr{" class=on" if start else ""}>' + "".join(
            f'<td{" class=r" if k in (0, 1, 2, 4, 5, 6, 8) else ""}>{esc(v)}</td>' for k, v in enumerate(tds))
            + "</tr>")
    return f'<table class="len">{head}<tbody>{"".join(body)}</tbody></table>'


def pool_chart(t, rows, L, lt, key):
    """The libraries of one pool, a row each: the P-sites per read length (start lengths dark, mixed hatched, the
    offset above each bar); with libraries.tsv, the regions and the frames of the start-length P-sites."""
    from matplotlib.figure import Figure
    n, side = len(rows), lt is not None
    h = 0.3 + 0.66 * n
    fig = Figure(figsize=(8.0, h))
    gs = fig.add_gridspec(n, 3 if side else 1, width_ratios=[1.9, 1, 1] if side else None, left=0.06, right=0.975,
                          top=1 - 0.16 / h, bottom=0.36 / h, hspace=0.85, wspace=0.08)
    L = np.array(L)
    tips = {}
    for i, (_, lib) in enumerate(rows.iterrows()):
        q = t[t["library"] == lib["library"]].set_index("length").reindex(L)
        ax = fig.add_subplot(gs[i, 0])
        top = max(float(q["share"].max()), 1e-9)
        for x, r in q.iterrows():
            if isna(r["share"]):
                continue
            start = bool(r["start_length"])
            kw = dict(width=0.78, linewidth=0, color=C["k-clean"] if start else C["k-noframe"])
            if r["mixed"]:
                kw.update(hatch="/////", edgecolor=C["surface"])
            b = ax.bar([x], [r["share"]], **kw)
            gid = f"{key}-{i}-{x}"
            b[0].set_gid(gid)
            tips[gid] = (f"{lib['name']}, {x} nt: {pct(r['share'])} of the P-sites"
                         + ("" if isna(r["offset"]) else f", offset {int(r['offset'])} nt")
                         + (f", mixed (nt {nt(r['beside_nt'])} has {f2(r['beside_share'])} of the peak)"
                            if r["mixed"] else "")
                         + (", start length" if start else f"; {r['qc']}"))
            if not isna(r["offset"]):
                ax.text(x, r["share"] + top * 0.06, f"{int(r['offset'])}", ha="center", va="bottom", fontsize=6,
                        color=C["ink2"] if start else C["muted"])
        ax.set_ylim(0, top * 1.45)
        ax.set_yticks([0, top], ["0", pct(top, 0)])
        ax.set_xlim(L[0] - 0.6, L[-1] + 0.6)
        ax.set_xticks(L if len(L) <= 16 else L[L % 2 == 0])
        if i < n - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("read length (nt)")
        ax.set_title(f"{lib['name']}  ·  {lib['kind']}", loc="left", fontsize=7, color=C["ink2"], pad=3)
        if side:
            r = lt.set_index("library").loc[lib["library"]]
            for col, parts, part, what, title in (
                    (1, REGION_VAR, 1, "", "start-length P-sites by region"),
                    (2, [(f"frame{f}", f"f{f}", f"frame {f}") for f in range(3)], 0, " inside the CDS",
                     "inside the CDS, by frame")):
                ax2 = fig.add_subplot(gs[i, col])
                share_bar(ax2, r, parts, part, what, lib["name"], f"{key}-{i}-s{col}", tips)
                if i == 0:
                    ax2.set_title(title, loc="left", fontsize=7, color=C["ink2"], pad=3)
                if i < n - 1:
                    ax2.set_xticks([])
                    ax2.spines["bottom"].set_visible(False)
                else:
                    ax2.set_xticks([0, 0.5, 1], ["0", "50%", "100%"])
    return svg(fig, key, f"P-sites per read length of the libraries of pool {rows['pool'].iat[0]}", tips)


def kernel_chart(k, pool, L, chosen, key):
    """Rows: the read lengths; cells: the kernel (harringtonine minus elongation at the annotated starts, a share of
    the row), binned."""
    from matplotlib.figure import Figure
    sub = k[(k["pool"] == pool) & k["length"].isin(L)]
    grid = sub.pivot(index="length", columns="nt", values="kernel").reindex(L).fillna(0)
    xs = grid.columns.to_numpy()
    h = 0.6 + 0.15 * len(L)
    fig = Figure(figsize=(3.9, h))
    fig.subplots_adjust(left=0.1, right=0.98, top=1 - 0.08 / h, bottom=0.42 / h)
    ax = fig.add_subplot()
    rows = np.arange(len(L))
    cells(ax, xs, rows, grid.to_numpy())
    ax.set_xlim(xs[0] - 0.5, xs[-1] + 0.5)
    ax.set_ylim(rows[-1] + 0.5, -0.5)
    ax.set_yticks(rows, [str(x) for x in L])
    for lab, x in zip(ax.get_yticklabels(), L):
        if x in chosen:
            lab.set_fontweight("bold")
            lab.set_color(C["ink"])
    ax.set_xticks([x for x in xs if x % 6 == 0], [nt(x) for x in xs if x % 6 == 0])
    ax.axvline(0, color=C["ink2"], linewidth=0.8, linestyle=(0, (2, 2)))
    ax.set_xlabel("P-site, nt from the start codon")
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    return svg(fig, key, f"Harringtonine minus elongation P-sites at the annotated starts per read length, pool {pool}")


def codon_chart(cod, classes, key):
    """Per class with calls: the share of each start codon among the tests and among the calls."""
    from matplotlib.figure import Figure
    fig = Figure(figsize=(8.0, 2.4))
    fig.subplots_adjust(left=0.06, right=0.99, top=0.88, bottom=0.22, wspace=0.2)
    tips = {}
    for j, (c, label) in enumerate(classes):
        ax = fig.add_subplot(1, len(classes), j + 1)
        sub = cod[cod["class"] == c].reset_index(drop=True)
        for i, r in sub.iterrows():
            for k, (what, colour, dx) in enumerate((("test_share", C["k-noframe"], -0.2),
                                                    ("call_share", C["k-clean"], 0.2))):
                b = ax.bar([i + dx], [r[what]], width=0.38, color=colour, linewidth=0)
                b[0].set_gid(f"{key}-{j}-{i}-{k}")
                tips[f"{key}-{j}-{i}-{k}"] = (f"{label}, {r['codon']}: {pct(r['call_share'])} of the calls "
                                              f"({count(r['called'])}), {pct(r['test_share'])} of the tests "
                                              f"({count(r['tests'])})")
        ax.set_xticks(range(len(sub)), sub["codon"], rotation=90, fontsize=6.5)
        ax.set_xlim(-0.6, len(sub) - 0.4)
        ax.set_ylim(0, max(float(sub[["test_share", "call_share"]].to_numpy().max()), 1e-9) * 1.1)
        ax.yaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%")
        ax.set_title(label, loc="left", fontsize=7.5, color=C["ink2"], pad=3)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    return svg(fig, key, "Start codons of the tests and of the calls per class", tips)


# ---------- text ----------

def start_command(stats, name):
    """The command of the run, from its inputs (the options that differ from the defaults)."""
    out = [f"pitstop starts --libraries {stats['libraries_sheet']} --gtf {stats['gtf']} --fasta {stats['fasta']}"]
    out.append("--scan" if stats["candidates_from"] == "scan" else f"--orfs {stats['candidates_from']}")
    if stats["start_lengths_given"]:
        out.append(f"--start-lengths {stats['start_lengths_given']}")
    if length_window(stats) != sx.WINDOW:
        out.append(f"--length-window {stats['length_window']}")
    if stats.get("start_codons") and tuple(stats["start_codons"].split(",")) != sx.START_CODONS:
        out.append(f"--start-codons {stats['start_codons']}")
    if float(stats["fdr"]) != sx.FDR:
        out.append(f"--fdr {stats['fdr']}")
    return " ".join(out + [f"--out-prefix {name}"])


def codon_classes(cod):
    """The classes of new starts (not the annotated ones) with calls: [(class, label)]."""
    if cod is None:
        return []
    return [(c, label) for c, label in CALL_CLASSES[1:] if (cod.loc[cod["class"] == c, "called"] > 0).any()]


def atg_line(cod, c, label):
    r = cod[(cod["class"] == c) & (cod["codon"] == "ATG")]
    if r.empty:
        return None
    r = r.iloc[0]
    return f"{label}: ATG is {pct(r['call_share'])} of the calls and {pct(r['test_share'])} of the tests"


def start_checks(d, chosen, mixed, paired, moved, pairs):
    """(name, state, value, basis, section) per check. Only the pairs have a cut of their own."""
    stats, t = d["stats"], d["t"]
    lo, hi = length_window(stats)
    passing = sorted(int(x) for x in t.loc[t["qc"] == "pass", "length"].unique())
    far = pairs[pairs["digestion"] == "warn"]
    weak = pairs[pairs["drug"] == "warn"]
    out = [("Same digestion", "warn" if len(far) else "pass",
            ", ".join(f"{r.label} ({dist(r.distance)})" for r in far.itertuples()) or f"every pair below "
                                                                                   f"{sx.PAIR_DISTANCE}",
            f"read-length distance of a pair < {sx.PAIR_DISTANCE}", "pairs")]
    if d["lt"] is not None:
        out.append(("Harringtonine effect", "warn" if len(weak) else "pass",
                    ("start peak not raised: " + ", ".join(weak["label"])) if len(weak) else
                    "start peak raised in every pair", "initiation start peak above its elongation library's",
                    "pairs"))
    out += [("Start lengths", "fail" if not chosen else "info" if stats["start_lengths_given"] else "pass",
             f"{span(chosen)} nt" if chosen else "none",
             f"given; the QC passes {span(passing)} nt" if stats["start_lengths_given"] else
             f"the QC, in the window {lo}–{hi} nt", "lengths"),
            ("Mixed lengths", "info" if mixed else "pass",
             ", ".join(f"{x} nt" + (" (pair)" if x in paired else "") for x in mixed) or "none",
             f"beside the start peak ≥ {sx.MIXED} of it, in a harringtonine library; a pair length in half of them "
             "or more", "lengths")]
    if chosen:
        moves = moved + moved_lengths(t, chosen, pair=True)
        out.append(("Calibration", "info" if moves else "pass",
                    "; ".join(f"{x} nt moved {nt(v)}" + (" (pair)" if x in paired else "") for x, v in sorted(moves))
                    or "no start length moved",
                    f"harringtonine start peak within {sx.MAX_SHIFT} nt of nt 0; a pair onto nt 0 and 1", "lengths"))
        key = f"q<{float(stats['fdr']):g} annotated"
        out.append(("Annotated starts", "info", f"{count(float(stats[key]))} of "
                    f"{count(float(stats['tests annotated']))} called", f"q < {stats['fdr']}", "calls"))
        classes = codon_classes(d["codons"])
        line = next((x for x in (atg_line(d["codons"], c, label) for c, label in classes) if x), None)
        if line:
            out.append(("Start codons", "info", line, "ATG, CTG and GTG start most; AAG, AGG and ATA least",
                        "calls"))
    return out


def start_steps(d, name, mixed, paired, moved, pairs):
    """(used, steps), each a list of (title, text) in order: what this run did, and what to do about it."""
    stats, t = d["stats"], d["t"]
    out, used = [], []
    lo, hi = length_window(stats)
    chosen = int_list(stats["start_lengths"])
    passing = sorted(int(x) for x in t.loc[t["qc"] == "pass", "length"].unique())
    far = pairs[pairs["digestion"] == "warn"]
    weak = pairs[pairs["drug"] == "warn"]
    if len(far):
        out.append((f"Check the pairs of another digestion: {', '.join(far['label'])}.",
                    f"Their two libraries differ in read lengths (distance {sx.PAIR_DISTANCE} or more): the initiation "
                    "library was digested differently from its elongation library, and the elongation reads then "
                    "expect another background. Pair the initiation library with an elongation library of the same "
                    "digestion, or leave the pair out."))
    if len(weak):
        out.append((f"Harringtonine did not raise the start peak: {', '.join(weak['label'])}.",
                    "In these pairs the initiation library has no higher start peak than its elongation library: the "
                    "drug did not stop the ribosomes at the starts. Leave the pair out, or make the library again."))
    if not chosen:
        out.append(("No read length passes the QC, so starts stopped.",
                    f"A start length must be in the window {lo}–{hi} nt, peak within {sx.MAX_SHIFT} nt of the "
                    "annotated starts in every harringtonine library, and have a clear frame in every library. Start "
                    "lengths says why each length fails. Choose the start lengths by hand with --start-lengths LO-HI, "
                    "or leave out the pair that fails, and run starts again."))
    elif stats["start_lengths_given"]:
        used.append((f"Start lengths given: {span(chosen)} nt.",
                     f"The QC passes {span(passing)} nt. Check the given lengths in each pair: the P-sites of a start "
                     "length peak at nt 0 of the annotated starts, in frame 0."))
    else:
        used.append((f"This run used the start lengths {span(chosen)} nt.",
                     f"Of the lengths in the window {lo}–{hi} nt, these are calibrated in every harringtonine "
                     "library, and have a clear frame in every library (the run of such lengths with the "
                     "most P-sites). Only their P-sites count. You do not need to pass them on; --start-lengths LO-HI "
                     "chooses them by hand."))
    if mixed:
        used.append((f"Mixed lengths: {span(mixed)} nt." + (" All are pair lengths." if paired == mixed else
                                                             f" Pair lengths: {span(paired)} nt." if paired else ""),
                     "In some harringtonine libraries the nucleotide beside the start peak holds half of the peak's "
                     "P-sites or more: two groups of 5′ ends, and one offset puts the smaller group 1 nt off. A length "
                     "that is mixed in half of the harringtonine libraries or more is a pair length: in every library, "
                     "starts anchors it on the peak and the nucleotide beside it, and moves the two onto nt 0 and 1 "
                     "of the start codon, the centre of a start's window. The inspect report of the library shows the "
                     "mixed lengths too."))
    if moved:
        used.append(("Some initiation libraries need other offsets.",
                     "The calibration moves the P-sites of " + ", ".join(f"{x} nt by {nt(v)}" for x, v in moved)
                     + " in some harringtonine libraries: their 5′ ends are longer or shorter than the elongation "
                     "libraries' offsets expect. starts corrects this per library and length; the tab of each pair "
                     "shows which."))
    if chosen and stats["candidates_from"] == "scan":
        out.append(("Count the ORFs of the catalogue in the elongation libraries.",
                    f"{name}.catalogue.tsv has the uORFs and uoORFs of the called starts, and every annotated CDS. "
                    "Count them with orfs --orfs on each elongation library, with the pool's offsets, then score "
                    "them: frame evidence comes from ribosomes that elongate."))
    elif chosen:
        out.append((f"The calls are in {name}.starts.tsv.",
                    f"A start with q < {stats['fdr']} has more initiation P-sites in its window than the elongation "
                    "libraries expect."))
    return used, out


# ---------- sections ----------

def more(summary, body):
    """Supplementary content, folded until opened."""
    return f'<details class="more"><summary>{esc(summary)}</summary>{body}</details>\n'


def pair_panel(d, r, L, chosen, charts, i):
    """The panel of one pair (the picker shows one at a time). charts: lengths, start (the start lengths, or every
    length of the window if there are none), other (the other lengths, or None), stop, cds."""
    t, lt = d["t"], d["lt"]
    items = [(r["digestion"], f"Read-length distance {dist(r['distance'])}: "
              + ("the two libraries have the same read lengths, as for one digestion." if r["digestion"] == "pass"
                 else f"{sx.PAIR_DISTANCE} or more, so the libraries were digested differently. Pair the initiation "
                      "library with an elongation library of the same digestion."))]
    if lt is not None:
        a, b = lt.set_index("library").loc[r["harringtonine"]], lt.set_index("library").loc[r["elongation"]]
        items.append((r["drug"],
                      "No start peak measured." if r["drug"] == "unknown" else
                      f"Start peak {fold(a['start_enrichment'])} in the initiation library against "
                      f"{fold(b['start_enrichment'])} in the elongation library: "
                      + ("harringtonine did not raise it. Leave the pair out, or make the library again."
                         if r["drug"] == "warn" else "harringtonine stopped the ribosomes at the starts.")))
        items.append(("info", f"The start lengths hold {pct(a['start_share'])} of the initiation library's P-sites "
                              f"and {pct(b['start_share'])} of the elongation library's."))
    moved = moved_lengths(t, chosen, r["harringtonine"])
    paired = moved_lengths(t, chosen, r["harringtonine"], pair=True)
    if paired:
        items.append(("info", "The calibration moves the pair lengths "
                      + ", ".join(f"{x} nt by {nt(v)}" for x, v in paired)
                      + " in the initiation library: their two groups of 5′ ends go onto nt 0 and 1."))
    if moved:
        items.append(("info", "The calibration moves the initiation library's "
                      + ", ".join(f"{x} nt by {nt(v)}" for x, v in moved)
                      + ": its 5′ ends differ from the elongation offsets."))
    h = t[(t["library"] == r["harringtonine"]) & t["length"].isin(L) & t["mixed"]]
    if len(h):
        items.append(("info", f"Mixed in the initiation library: {span(sorted(int(x) for x in h['length']))} nt."))
    both = names(t[t["library"].isin([r["harringtonine"], r["elongation"]]) & t["length"].isin(L)], "library")
    frames = legend([(None, f"var(--f{f})", f"frame {f}") for f in range(3)])
    profile_note = ("Initiation above, elongation below (faded). The two libraries have different y scales: the tick "
                    "of each gives its highest bar, so compare the shapes, not the heights. P-sites from the pool's "
                    "offsets, before the calibration: a length that the calibration moves peaks off nt 0. Colour: "
                    "the frame of the CDS.")
    other = [x for x in L if x not in chosen]
    return (f'<div class="libpanel" id="pair-{i}" role="region" aria-label="{esc(r["label"])}"'
            f'{" hidden" if i else ""}>'
            f'<h3>{esc(r["label"])}</h3><p class="lead">Initiation: <code>{esc(r["h_name"])}</code>. Elongation: '
            f'<code>{esc(r["e_name"])}</code>.</p>'
            + findings(items)
            + fig_block("Read lengths and frames of the two libraries", charts["lengths"],
                        f"Above: the share of each library's P-sites inside the CDS ({sx.INTERIOR} nt or more from "
                        "both ends) per read length: the footprints of elongating ribosomes, so the lengths that the "
                        "digestion leaves. The distance is half the summed difference of the two. Below: the frames "
                        "of those P-sites per read length, initiation left, elongation right (faded). Where the offset "
                        "is right, frame 0 is the largest.",
                        legend([(None, "var(--k-clean)", "initiation"), (None, "var(--k-short)", "elongation")]
                               + [(None, f"var(--f{f})", f"frame {f}") for f in range(3)]))
            + fig_block(f"P-sites at the annotated starts, the start lengths ({span(chosen)} nt)" if chosen else
                        "P-sites at the annotated starts, per read length", charts["start"],
                        "A row per read length. " + profile_note, frames)
            + (more(f"The other read lengths of the window: {span(other)} nt",
                    fig_block("P-sites at the annotated starts, the other read lengths", charts["other"],
                              "As above. A length is mixed where the nt beside its start peak holds half of the "
                              "peak's P-sites or more.", frames)) if chosen and other else "")
            + more("The pair per read length, in numbers",
                   '<div class="scroll">' + length_table(t, r, L)
                   + '<p class="note">Per read length in the window. P-sites: the share of the library\'s P-sites. '
                     "Offset: nt from the 5′ end to the P-site, after the calibration (in brackets, its move). Start "
                     "peak: the nt of the most P-sites at the annotated starts, before the calibration; mixed: the nt "
                     "beside it holds this share of the peak or more. Top phase: the share of the top frame inside "
                     "the CDS. Start lengths in bold.</p></div>\n" + tsv_view(both, f"{d['name']}.lengths.tsv"))
            + more("The stop codons and along the CDS, per read length",
                   fig_block("P-sites at the stop codons, per read length", charts["stop"],
                             "As at the starts: initiation above, elongation below, different y scales. The P-site "
                             "of a terminating ribosome is the last sense codon, nt −3: a peak there in the "
                             "elongation library confirms its offset.", frames)
                   + fig_block("Along the CDS, per read length", charts["cds"],
                               f"The P-sites in {sx.CDS_BINS} bins from the start to the stop codon, as a share of "
                               "the library's P-sites in the CDS. In an initiation library the first bin holds the "
                               "start peak, and after it the line climbs towards the stop: the ribosomes ran off the "
                               "5′ part first.",
                               legend([(None, "var(--ink)", "initiation"), ("line", None, "elongation")])))
            + "</div>\n")


def calls_section(d, name):
    stats, cod = d["stats"], d["codons"]
    fdr = stats["fdr"]
    key = f"q<{float(fdr):g}"
    num = lambda k: float(stats[k]) if stats.get(k, "") != "" else np.nan
    scan = stats["candidates_from"] == "scan"
    rows = [(label, count(num(f"tests {c}")), count(num(f"candidates {c}")), count(num(f"{key} {c}")),
             count(num(f"expected_false {c}"))) for c, label in CALL_CLASSES if f"tests {c}" in stats]
    table = pd.DataFrame(rows, columns=["Class", "Tests", "Candidates", f"q < {fdr}", "Expected false"])
    tl = [("Candidates", count(num("candidates")) if scan else count(num("orf_table_rows")),
           "start codons with initiation reads" if scan else "rows of the ORF table"),
          ("Annotated starts", count(num(f"{key} annotated")), f"of {count(num('tests annotated'))}, q < {fdr}"),
          ("Null windows", count(num("null_windows")), "no start codon nearby")]
    if scan:
        tl.insert(1, ("Called", count(num("called")), f"q < {fdr}, not stopped"))
    stopped = (f" {count(num('stopped'))} candidates with q < {fdr} were not called: a called or annotated start stops "
               f"the candidates within {sx.NEAR[1]} nt and in its codons 1–10, where harringtonine raises the reads."
               if scan else "")
    body = (tiles(tl)
            + f'<div class="panel scroll">{html_table(table, ["Tests", "Candidates", f"q < {fdr}", "Expected false"])}'
            '<p class="note">A class is one BH family: the annotated starts, and the other candidates by region. Tests: '
            "the start codons of the class (a scan tests those with harringtonine reads; the others count with p = 1). "
            "Expected false: the share of the null windows (windows with no start codon nearby) at the p of the last "
            "call, times the tests; where they are a large part of the calls of a class, use those calls with care."
            f"{stopped}</p></div>\n")
    if scan:
        orfs = pd.DataFrame([(ty, count(num(f"start_orfs {ty}")), count(num(f"catalogue {ty}")))
                             for ty in ORF_TYPES if f"start_orfs {ty}" in stats or f"catalogue {ty}" in stats],
                            columns=["Type", "start_orfs.tsv", "catalogue.tsv"])
        body += (f'<div class="panel scroll narrow">{html_table(orfs, ["start_orfs.tsv", "catalogue.tsv"])}'
                 '<p class="note">The ORFs of the called starts, to the first in-frame stop, by type against the '
                 f"annotated CDS. The catalogue ({esc(name)}.catalogue.tsv): their uORFs and uoORFs, and every "
                 "annotated CDS.</p></div>\n")
    classes = codon_classes(cod)
    if classes:
        lines = [x for x in (atg_line(cod, c, label) for c, label in classes) if x]
        body += (findings([("info", "Real starts are rich in ATG, then CTG and GTG; AAG, AGG and ATA start least "
                                    "often. Calls whose start codons look like those of the tests are mostly noise: "
                                    "restrict --start-codons, or lower --fdr.")]
                          + [("info", f"{x}.") for x in lines])
                 + fig_block("Start codons of the calls against the tests", d["codon_chart"],
                             "Per class of new starts: the share of each start codon among the tests (all such codons "
                             "of the class) and among the calls.",
                             legend([(None, "var(--k-noframe)", "tests"), (None, "var(--k-clean)", "calls")]),
                             tsv_view(cod, f"{name}.start_codons.tsv")))
    return body


def pool_totals(stats):
    """The P-sites of the start lengths per pool, as a panel."""
    num = lambda k: float(stats[k]) if stats.get(k, "") != "" else np.nan
    pools = sorted({k.split(" ", 1)[1] for k in stats if k.startswith("harringtonine_psites ")})
    totals = pd.DataFrame([(pl, count(num(f"harringtonine_psites {pl}")), count(num(f"elongation_psites {pl}")))
                           for pl in pools], columns=["Pool", "Harringtonine P-sites", "Elongation P-sites"])
    return (f'<div class="panel scroll narrow">{html_table(totals, ["Harringtonine P-sites", "Elongation P-sites"])}'
            '<p class="note">The P-sites of the start lengths per pool, after the calibration: the counts of the '
            "test.</p></div>\n")


def starts_report(prefix):
    """<prefix>.starts.html from the tables of a starts run."""
    import matplotlib
    d = load_starts(prefix)
    stats, t, lt, kern = d["stats"], d["t"], d["lt"], d["kern"]
    name = d["name"] = os.path.basename(prefix)
    libs = library_order(t)
    pairs = pair_order(d["pairs"], lt)
    L = window_lengths(t, stats)
    chosen = int_list(stats["start_lengths"])
    lo, hi = length_window(stats)
    harr = t[(t["kind"] == "harringtonine") & t["length"].isin(L)]
    mixed = sorted(int(x) for x in harr.loc[harr["mixed"], "length"].unique())
    paired = sorted(int(x) for x in harr.loc[harr["pair"], "length"].unique())
    moved = moved_lengths(t, chosen)
    pools = list(dict.fromkeys(libs["pool"]))
    first = [x for x in L if x in chosen] or L                   # the start lengths show, the others are folded
    other = [x for x in L if x not in first]
    with matplotlib.rc_context(RC):
        overview = pairs_chart(pairs, lt, "pairs")
        pair_charts = [dict(lengths=lengths_pair_chart(t, r, f"plen{i}"),
                            start=profile_rows(t, d["prof"], r, first, "start", f"pstart{i}"),
                            other=profile_rows(t, d["prof"], r, other, "start", f"pother{i}") if other else None,
                            stop=profile_rows(t, d["prof"], r, L, "stop", f"pstop{i}"),
                            cds=cds_grid(t, d["prof"], r, L, f"pcds{i}"))
                       for i, r in pairs.iterrows()]
        pool_charts = {pl: pool_chart(t, libs[libs["pool"] == pl], L, lt, f"pool{i}") for i, pl in enumerate(pools)}
        kernel_charts = ({pl: kernel_chart(kern, pl, L, chosen, f"kern{i}") for i, pl in enumerate(pools)}
                         if kern is not None else {})
        classes = codon_classes(d["codons"])
        d["codon_chart"] = codon_chart(d["codons"], classes, "codons") if classes else ""

    files = [f"{name}.{x}" for x in ("starts_stats.tsv", "pairs.tsv", "lengths.tsv", "profiles.tsv")]
    if chosen:
        files += [f"{name}.{x}" for x in ("libraries.tsv", "kernel.tsv", "start_codons.tsv", "starts.tsv")]
        if stats["candidates_from"] == "scan":
            files += [f"{name}.{x}" for x in ("start_orfs.tsv", "catalogue.tsv")]
    checks = start_checks(d, chosen, mixed, paired, moved, pairs)
    worst = lambda r: next(s for s in ("fail", "warn", "info", "unknown", "pass") if s in (r["digestion"], r["drug"]))
    side = (f"<div>{LOCKUP}</div>\n<div><div class=\"kind\">starts report</div><div class=\"lib\">{esc(name)}</div>"
            "</div>\n<div><h3>Checks</h3><ul class=\"flaglist\">"
            + "".join(f'<li><a href="#{sec}">{icon(state)}<span>{esc(label)}</span></a></li>'
                      for label, state, _, _, sec in checks)
            + "</ul></div>\n<div><h3>Pairs</h3><ul class=\"flaglist\">"
            + "".join(f'<li><a href="#pair" data-pair="pair-{i}">{icon(worst(r))}<span>{esc(r["label"])}</span></a></li>'
                      for i, r in pairs.iterrows())
            + "</ul></div>\n<div><h3>Sections</h3><ul class=\"toc\">"
            + "".join(f'<li><a href="#{sid}">{esc(title)}</a></li>' for sid, title in S_SECTIONS)
            + "</ul></div>\n<div class=\"files\"><h3>Tables</h3>"
            + "".join(f"<code>{esc(f)}</code>" for f in files) + "</div>\n")

    n_pools = stats["pools"]
    head = ('<header><div class="eyebrow">pitstop starts · initiation against elongation</div>'
            f"<h1>{esc(name)}</h1><dl class=\"meta\">"
            f"<dt>Command</dt><dd><code>{esc(start_command(stats, name))}</code></dd>"
            f"<dt>pitstop</dt><dd>{esc(stats['pitstop'])}</dd>"
            f"<dt>Libraries</dt><dd>{esc(stats['libraries'])} in {esc(n_pools)} "
            f"pool{'' if n_pools == '1' else 's'}; {count(float(stats['transcripts_with_cds']))} of "
            f"{count(float(stats['transcripts']))} transcripts with an annotated CDS</dd>"
            f"<dt>Start lengths</dt><dd>"
            + (f"{span(chosen)} nt" + (" (given)" if stats["start_lengths_given"] else "") if chosen else "none")
            + f" (window {lo}–{hi} nt)</dd></dl></header>\n")

    used, todo = start_steps(d, name, mixed, paired, moved, pairs)
    items = lambda xs: "".join(f'<li><div class="t">{esc(title)}</div><p>{esc(text)}</p></li>' for title, text in xs)
    steps = (f'<h3>Used in this run</h3><ul class="steps">{items(used)}</ul>' if used else "") \
        + f'<h3>Next steps</h3><ol class="steps">{items(todo)}</ol>'
    check_html = ('<table class="flags"><thead><tr><th>Check</th><th>State</th><th>Value</th><th>Basis</th><th></th>'
                  "</tr></thead><tbody>" + "".join(
                      f"<tr><td>{esc(label)}</td><td>{state_cell(state)}</td><td>{esc(value)}</td>"
                      f'<td>{esc(basis)}</td><td><a href="#{sec}">details</a></td></tr>'
                      for label, state, value, basis, sec in checks) + "</tbody></table>")
    verdict = section("verdict", "Verdict and next steps", (
        f'<div class="vgrid"><div class="panel nextbox">{steps}</div>'
        f'<div class="panel scroll">{check_html}<p class="note">The pairs have a cut of their own; the QC of the read '
        "lengths decides the start lengths, and the null windows estimate the false calls. A warning does not change "
        "the exit code.</p></div></div>\n"), fold=False)

    far, weak = (pairs["digestion"] == "warn").any(), (pairs["drug"] == "warn").any()
    pair_find = [("warn" if far else "pass",
                  "Some pairs differ in read lengths: another digestion." if far else
                  f"Every pair has a read-length distance below {sx.PAIR_DISTANCE}: one digestion.")]
    if lt is not None:
        pair_find.append(("warn" if weak else "pass",
                          "In some pairs harringtonine did not raise the start peak." if weak else
                          "In every pair the initiation library has the higher start peak."))
    pairs_sec = section("pairs", "Pairs", (
        '<p class="lead">A pair is a row of the library sheet: an initiation library and the elongation library that '
        "gives its background. The two must have the same digestion, and the drug must raise the start peak of the "
        'initiation library. Each pair in detail: <a href="#pair">Pair detail</a>.</p>\n'
        + findings(pair_find)
        + fig_block("The pairs at a glance", overview,
                    f"Left: the read-length distance of the CDS-interior P-sites (0: the same lengths, 1: none in "
                    f"common); dashed: the cut, {sx.PAIR_DISTANCE}. Right: the start peak, the P-sites at the "
                    "annotated start codons per start over the P-sites per in-frame codon inside the CDS, for the "
                    "start lengths after the calibration. Harringtonine raises it.",
                    legend([(None, "var(--k-short)", "distance"), (None, "var(--warn)", "distance ≥ cut: warn")]
                           + ([(None, "var(--k-clean)", "initiation"), (None, "var(--k-short)", "elongation")]
                              if lt is not None else [])),
                    tsv_view(names(d["pairs"], "harringtonine", "elongation"), f"{name}.pairs.tsv"))))

    per_length = t[t["length"].isin(L)].drop_duplicates("length")[["length", "start_length", "qc"]]
    on = [i for i, x in enumerate(per_length["start_length"]) if x]
    per_length = per_length.assign(length=[f"{x} nt" for x in per_length["length"]],
                                   start_length=per_length["start_length"].map({True: "yes", False: ""}))
    per_length.columns = ["Length", "Start length", "QC"]
    state_len = next(c[1] for c in checks if c[0] == "Start lengths")
    len_find = [(state_len, f"Start lengths: {span(chosen)} nt." if chosen else "No read length passes the QC.")]
    leg = legend([(None, "var(--k-clean)", "start length"), (None, "var(--k-noframe)", "other length")]
                 + [(None, "repeating-linear-gradient(135deg, var(--k-noframe) 0 2px, var(--surface) 2px 3.5px)",
                     "mixed")]
                 + ([(None, f"var(--{var})", label) for _, var, label in REGION_VAR]
                    + [(None, f"var(--f{f})", f"frame {f}") for f in range(3)] if lt is not None else []))
    blocks = "".join(fig_block(f"Pool {pl}", pool_charts[pl],
                               "Left: the P-sites per read length, as a share of the library's; above each bar its "
                               "offset (nt from the 5′ end to the P-site, after the calibration)."
                               + (" Middle and right: the P-sites of the start lengths, after the calibration, by "
                                  "region, and by frame inside the CDS." if lt is not None else ""),
                               leg, tsv_view(names(lt[lt["pool"] == pl], "library"), f"{name}.libraries.tsv")
                               if lt is not None else "")
                     for pl in pools)
    if kern is not None:
        kblocks = "".join(fig_block(f"Pool {pl}", kernel_charts[pl], "", small=True) for pl in pools)
        kernel = more("Per pool: harringtonine at the annotated starts", (
            '<p class="lead">Per pool and read length: the harringtonine P-sites minus the elongation P-sites (scaled '
            "to match on codons 11–20), as a share of nt −15 to +32 from the annotated starts, after the calibration. "
            "Where harringtonine works and the offset is right, most of a row is at nt 0. The start lengths are in "
            "bold.</p>\n"
            + legend([(None, f"var(--h{i})", f"{'< ' + str(HEAT_BINS[0]) if i == 0 else '≥ ' + str(HEAT_BINS[i - 1])}")
                      for i in range(8)])
            + f'<div class="figrow fill">{kblocks}</div>\n'
            + f'<p class="note">The table: <code>{esc(name)}.kernel.tsv</code>.</p>\n'))
        calls = section("calls", "Calls", calls_section(d, name))
    else:
        kernel = ('<p class="note">Harringtonine at the annotated starts per pool and read length: not computed, no '
                  "start lengths.</p>\n")
        calls = section("calls", "Calls", '<p class="lead">No calls: no start lengths.</p>\n')
    lengths = section("lengths", "Start lengths", (
        findings(len_find)
        + f'<div class="panel scroll">{html_table(per_length, ["Length"], cls="len", on=on)}<p class="note">The QC of '
          f"each read length in the window, over all libraries. A start length is in the window {lo}–{hi} nt, its "
          f"P-sites peak within {sx.MAX_SHIFT} nt of the annotated starts ({sx.MIN_ANCHOR} P-sites or more there; "
          "for a pair length, the nt beside the peak too) in every harringtonine library, and "
          f"{sx.MIN_FRAME} P-sites or more inside the CDS have {sx.MIN_TOP} or more in their top frame in every "
          "library. Of the runs of such lengths, the one with the most P-sites.</p></div>\n"
        + more("Per pool: the read lengths, offsets, regions and frames of each library",
               blocks + (pool_totals(stats) if chosen else ""))
        + kernel))

    options = "".join(
        f'<optgroup label="{esc(pl)}">'
        + "".join(f'<option value="pair-{i}"{" selected" if i == 0 else ""}>{esc(r["label"])}'
                  + ("" if worst(r) == "pass" else f" ({STATE[worst(r)][1].lower()})") + "</option>"
                  for i, r in pairs.iterrows() if r["pool"] == pl)
        + "</optgroup>" for pl in dict.fromkeys(pairs["pool"]))
    picker = ('<div class="pairbar"><label for="pair-pick">Pair</label>'
              '<button type="button" data-step="-1" aria-label="Previous pair">‹</button>'
              f'<select id="pair-pick" autocomplete="off">{options}</select>'
              '<button type="button" data-step="1" aria-label="Next pair">›</button>'
              "<span>all of this section is for the chosen pair</span></div>\n")
    panels = "".join(pair_panel(d, r, L, chosen, pair_charts[i], i) for i, r in pairs.iterrows())
    detail = section("pair", "Pair detail", (
        '<p class="lead">One pair at a time: its two libraries per read length, at the annotated starts and along '
        "the CDS. The sections above are for the whole run.</p>\n" + picker + panels))

    scope = section("scope", "Not in this report", (
        '<ul class="scope">'
        "<li><strong>The QC of each library:</strong> read ends, untemplated nucleotides, under-digestion, "
        "complexity. Run <code>pitstop inspect</code> on each library.</li>"
        "<li><strong>The agreement of the replicates:</strong> starts sums the pairs of a pool for the test.</li>"
        "<li><strong>The frame evidence of the new ORFs:</strong> <code>pitstop orfs</code> and "
        "<code>pitstop score</code> on the elongation libraries.</li>"
        "</ul>\n"), fold=False)
    foot = ("<footer><div>Tables: " + ", ".join(f"<code>{esc(f)}</code>" for f in files) + "</div><div>Each number "
            "in this report is also in a table. A rerun with the same inputs gives a byte-identical report.</div>"
            "</footer>\n")
    out = page(f"{name} · pitstop starts", side, head + verdict + calls + pairs_sec + lengths + detail + scope + foot,
               "Checks and contents")
    path = f"{prefix}.starts.html"
    with open(path, "w", encoding="utf-8") as f:
        f.write(out)
    log.info("wrote %s", path)
    return path

"""QC figures of a `pitstop starts` run, from <prefix>.lengths.tsv and <prefix>.profiles.tsv: which read lengths each
library has, where their P-sites fall, and which lengths the QC takes as start lengths. Needs matplotlib.

Page 1, a row per library: its P-sites per read length (the start lengths dark, the window shaded), and the regions
and CDS frames of its start-length P-sites (the frames after the calibration). Then a page per library, a row per read
length: the frames inside the CDS, the P-sites at the annotated starts and at the stop codons (per nt, coloured by
frame), and along the CDS (bins, per frame). These P-sites are from the pool's offsets, before the calibration: a
length that the calibration moves peaks off nt 0, and its row says by how much."""
import logging
import os
import textwrap

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

FRAME = ("#2a78d6", "#eb6834", "#1baf7a")            # frame 0, 1, 2
REGION = {"leader": "#eda100", "CDS": "#008300", "trailer": "#e87ba4", "transcript": "#6250d6"}
REGION_LABEL = {"leader": "leader", "CDS": "CDS", "trailer": "trailer", "transcript": "no CDS"}
INK, INK2, MUTED, AXIS, SHADE = "#0b0b0b", "#52514e", "#898781", "#c3c2b7", "#efeee9"
A4 = (8.27, 11.69)
RC = {"font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
      "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
      "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2, "ytick.major.size": 2,
      "axes.labelcolor": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
      "pdf.fonttype": 42, "font.family": "DejaVu Sans"}


def compact(v):
    """A tick label: 950, 12k, 1.2M."""
    return f"{v / 1e6:.1f}M" if v >= 1e6 else f"{v / 1e3:.0f}k" if v >= 1e4 else f"{v / 1e3:.1f}k" if v >= 1e3 \
        else f"{v:.0f}"


def spans(lengths):
    """24,25,26,28 -> '24-26, 28'."""
    out, run = [], []
    for x in lengths:
        if run and x != run[-1] + 1:
            out.append(run)
            run = []
        run.append(x)
    out += [run] if run else []
    return ", ".join(f"{r[0]}-{r[-1]}" if len(r) > 1 else f"{r[0]}" for r in out) or "none"


def short_reasons(qc):
    """lengths.tsv's qc column, one short line per reason, for a row label."""
    out = []
    for r in qc.split("; "):
        w = r.split()
        if r.startswith("outside the window"):
            out.append(f"outside {w[-1]}")
        elif r.startswith("outside the run"):
            out.append("outside the run")
        elif r.startswith("not calibrated"):
            out.append(f"not calibrated in {w[3]}/{w[5]} harr.")
        elif r.startswith("no clear frame"):
            out.append(f"no clear frame in {w[4]}/{w[6]}")
        else:
            out.append(textwrap.shorten(r, 30))
    return out


def signed(v):
    return f"{v:+d}" if v else "0"


def libraries(qc):
    """One row per library: pool, then harringtonine before elongation, then the order of the table."""
    libs = qc.drop_duplicates("library")[["library", "pool", "kind"]].reset_index(drop=True)
    libs["order"] = libs.index
    libs["elong"] = libs["kind"] != "harringtonine"
    return libs.sort_values(["pool", "elong", "order"], kind="stable").reset_index(drop=True)


def selected(qc):
    """The start lengths; without them the lengths in the window."""
    chosen = qc.loc[qc["start_length"], "length"].unique()
    if len(chosen):
        return sorted(int(x) for x in chosen), "start lengths"
    inside = ~qc["qc"].str.startswith("outside the window")
    return sorted(int(x) for x in qc.loc[inside, "length"].unique()), "no start lengths: the lengths in the window"


def stacked(ax, shares, colours, labels, label_part):
    """Horizontal stacked bars, one per row of shares (rows x parts, each row sums to 1), top row first; the share
    of part label_part written inside its segment."""
    n = len(shares)
    left = np.zeros(n)
    for j, (c, lab) in enumerate(zip(colours, labels)):
        ax.barh(np.arange(n), shares[:, j], left=left, height=0.62, color=c, edgecolor="white", linewidth=0.8,
                label=lab)
        if j == label_part:
            for i in range(n):
                if shares[i, j] >= 0.12:
                    ax.text(left[i] + shares[i, j] / 2, i, f"{shares[i, j]:.0%}", ha="center", va="center",
                            color="white", fontsize=5.5, fontweight="bold")
        left += shares[:, j]
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1], ["0", "50%", "100%"])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(labels), frameon=False, fontsize=5.5,
              handlelength=1, handleheight=0.8, columnspacing=0.8, borderaxespad=0.3)


def overview(fig, qc, libs, name):
    """Page 1: the read lengths of every library, and the regions and frames of the start-length P-sites."""
    n = len(libs)
    lengths = np.sort(qc["length"].unique())
    chosen, what = selected(qc)
    inside = sorted(qc.loc[~qc["qc"].str.startswith("outside the window"), "length"].unique())
    hspace = 0.55
    gs = fig.add_gridspec(n, 3, left=0.05, right=0.98, top=0.9, bottom=0.05, width_ratios=[1.7, 1, 1],
                          hspace=hspace, wspace=0.12)
    fig.text(0.05, 0.965, f"Read lengths of {name}", fontsize=10, fontweight="bold", color=INK)
    fig.text(0.05, 0.945, f"{what}: {spans(chosen)}. Left: P-sites per read length, as a share of the library's; "
                          "the window shaded, the start lengths dark.", fontsize=6.5, color=INK2)
    fig.text(0.05, 0.932, "Middle and right: the start-length P-sites by region, and by frame inside the CDSs "
                          "(after the calibration). Tables: <prefix>.lengths.tsv.", fontsize=6.5, color=INK2)
    regions, frames = [], []
    for i, lib in libs.iterrows():
        q = qc[qc["library"] == lib["library"]].set_index("length").reindex(lengths)
        ax = fig.add_subplot(gs[i, 0])
        if inside:
            ax.axvspan(inside[0] - 0.5, inside[-1] + 0.5, color=SHADE, zorder=0, linewidth=0)
        ax.bar(lengths, q["share"].fillna(0), width=0.8, zorder=2,
               color=[INK2 if x in chosen else AXIS for x in lengths])
        top = max(float(q["share"].max()), 1e-9)
        ax.set_ylim(0, top * 1.15)
        ax.set_yticks([0, top], ["0", f"{top:.0%}"])
        ax.set_xlim(lengths[0] - 0.6, lengths[-1] + 0.6)
        ax.set_xticks(lengths[lengths % 2 == 0])
        if i < n - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("read length (nt)")
        ax.set_title(f"{os.path.basename(lib['library'])}  ·  {lib['kind']}, pool {lib['pool']}", loc="left",
                     fontsize=5.5, color=INK2, pad=2)
        s = q.loc[chosen]
        r = s[list(REGION)].sum().to_numpy(dtype=float)
        regions.append(r / max(r.sum(), 1))
        f = np.zeros(3)
        for length, row in s.iterrows():   # the calibration moves P-sites by `shift`: frame f -> f + shift
            f += np.roll(row[["frame0", "frame1", "frame2"]].to_numpy(dtype=float), int(row["shift"]))
        frames.append(f / max(f.sum(), 1))
    a = 1 / (2 * (1 + hspace))   # the spanning axes put row i's bar at the centre of row i's small multiple
    for col, shares, colours, labels, part in (
            (1, regions, list(REGION.values()), [REGION_LABEL[r] for r in REGION], 1),
            (2, frames, FRAME, ["frame 0", "frame 1", "frame 2"], 0)):
        ax = fig.add_subplot(gs[:, col])
        stacked(ax, np.array(shares), colours, labels, part)
        ax.set_ylim(n - 1 + a, -a)
        ax.set_xlabel("start-length P-sites by region" if col == 1 else "start-length P-sites in the CDS, by frame")


def library_page(fig, qc, prof, lib):
    """A page per library: per read length, the frames, the P-sites at the starts and stops, along the CDS."""
    q = qc[qc["library"] == lib["library"]].set_index("length").sort_index()
    p = prof[prof["library"] == lib["library"]]
    tables = {a: g.pivot(index="length", columns="x", values="psites").reindex(q.index, fill_value=0)
              for a, g in p[p["anchor"] != "cds"].groupby("anchor")}
    cds = p[p["anchor"] == "cds"]
    n = len(q)
    gs = fig.add_gridspec(n, 5, left=0.02, right=0.985, top=0.9, bottom=0.045,
                          width_ratios=[1.15, 0.5, 2.4, 1.9, 1.15], hspace=0.5, wspace=0.2)
    fig.text(0.02, 0.965, os.path.basename(lib["library"]), fontsize=10, fontweight="bold", color=INK)
    fig.text(0.02, 0.945, f"{lib['kind']}, pool {lib['pool']}. A row per read length, the start lengths in bold.",
             fontsize=6.5, color=INK2)
    fig.text(0.02, 0.932, "P-sites from the pool's offsets, before the calibration; colour: the frame of the CDS.",
             fontsize=6.5, color=INK2)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=c, label=f"frame {i}") for i, c in enumerate(FRAME)], loc="upper right",
               bbox_to_anchor=(0.985, 0.975), ncol=3, frameon=False, fontsize=6, handlelength=1, handleheight=0.8)
    for i, (length, row) in enumerate(q.iterrows()):
        last = i == n - 1
        lab = fig.add_subplot(gs[i, 0])
        lab.axis("off")
        start = bool(row["start_length"])
        lab.text(0, 0.95, f"{length} nt", fontsize=7.5, fontweight="bold" if start else "normal",
                 color=INK if start else MUTED, va="top", transform=lab.transAxes)
        note = [f"{row['share']:.1%} of P-sites"]
        if row["shift"]:
            note.append(f"calibration: {signed(int(row['shift']))} nt")
        note += ["start length"] if start else short_reasons(str(row["qc"]))
        lab.text(0, 0.6, "\n".join(note), fontsize=5, color=INK2 if start else MUTED, va="top", linespacing=1.1,
                 transform=lab.transAxes)

        ax = fig.add_subplot(gs[i, 1])
        fr = row[["frame0", "frame1", "frame2"]].to_numpy(dtype=float)
        fr = fr / max(fr.sum(), 1)
        ax.bar(range(3), fr, color=FRAME, width=0.75)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1], ["0", "1"])
        ax.set_xticks([])
        ax.text(2.6, 1.0, f"{fr.max():.2f}", ha="right", va="top", fontsize=5, color=INK2)
        if i == 0:
            ax.set_title("frame, CDS", loc="left")

        for col, anchor, title, xlabel in ((2, "start", "P-sites at the annotated starts", "nt from the start codon"),
                                           (3, "stop", "P-sites at the stop codons", "nt from the stop codon")):
            ax = fig.add_subplot(gs[i, col])
            t = tables[anchor]
            x = t.columns.to_numpy()
            y = t.loc[length].to_numpy(dtype=float)
            ax.axvline(0, color=AXIS, linewidth=0.6, linestyle=(0, (2, 2)), zorder=0)
            ax.bar(x, y, width=0.85, color=[FRAME[v % 3] for v in x], linewidth=0)
            top = max(y.max(), 1)
            ax.set_ylim(0, top * 1.3)
            ax.set_yticks([0, top], ["0", compact(top)])
            ax.set_xlim(x[0] - 1, x[-1] + 1)
            ax.set_xticks([v for v in x if v % 15 == 0])
            if anchor == "start" and row["at_peak"]:
                ax.text(x[-1], top * 1.3, f"peak at {signed(int(row['start_peak']))}", ha="right", va="top",
                        fontsize=5, color=INK2)
            if not last:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(xlabel)
            if i == 0:
                ax.set_title(title, loc="left")

        ax = fig.add_subplot(gs[i, 4])
        c = cds[cds["length"] == length].pivot(index="x", columns="frame", values="psites").reindex(columns=range(3),
                                                                                                    fill_value=0)
        share = c.to_numpy(dtype=float) / max(c.to_numpy().sum(), 1)
        for f in range(3):
            ax.plot(c.index, share[:, f], color=FRAME[f], linewidth=1)
        top = max(share.max(), 1e-9)
        ax.set_ylim(0, top * 1.15)
        ax.set_yticks([0, top], ["0", f"{top:.0%}"])
        ax.set_xticks([])
        if last:
            ax.set_xlabel("along the CDS, 5' to 3'")
        if i == 0:
            ax.set_title("P-sites along the CDS", loc="left")


def run(prefix, out=None):
    """<prefix>.lengths.pdf (or `out`) from <prefix>.lengths.tsv and <prefix>.profiles.tsv of `pitstop starts`."""
    try:
        import matplotlib
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.figure import Figure
    except ImportError:
        raise ImportError("pitstop plot needs matplotlib: pip install 'pitstop[plot]'") from None
    logging.getLogger("fontTools").setLevel(logging.WARNING)   # it logs every font subset at INFO
    out = out or f"{prefix}.lengths.pdf"
    qc = pd.read_csv(f"{prefix}.lengths.tsv", sep="\t", dtype={"library": str, "pool": str, "qc": str})
    prof = pd.read_csv(f"{prefix}.profiles.tsv", sep="\t", dtype={"library": str})
    libs = libraries(qc)
    with matplotlib.rc_context(RC), PdfPages(out, metadata={"Creator": "pitstop", "CreationDate": None}) as pdf:
        fig = Figure(figsize=A4)
        overview(fig, qc, libs, os.path.basename(prefix))
        pdf.savefig(fig)
        for _, lib in libs.iterrows():
            fig = Figure(figsize=A4)
            library_page(fig, qc, prof, lib)
            pdf.savefig(fig)
    log.info("wrote %s (%d pages)", out, 1 + len(libs))
    return out

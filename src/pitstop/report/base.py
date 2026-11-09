"""The frame of the HTML reports: the theme, the page, the figures as inline SVG, numbers as text, and the parts of a
page. Each report has its module (inspect.py, starts.py); page.html, style.css and script.js are files of this package.

A page needs no other file and no network. The figures are matplotlib SVG in the page. They are drawn in the light
colours of THEME, and then each of these colours becomes its CSS variable, so the figures follow the light or dark
theme of the browser with the rest of the page. No date and no path in a page: a rerun gives a byte-identical file.
Every number on a page comes from a table of the run; a page formats numbers, it does not calculate them.
"""
import html
import io
import logging
import re
from importlib.resources import files
from string import Template

import numpy as np
import pandas as pd

from .. import inspect as qc

log = logging.getLogger(__name__)

# CSS variable: (light, dark). FIGURE: the variables a figure can use (their light values must be unique).
THEME = {
    "bg": ("#f9f9f7", "#0d0d0d"), "surface": ("#fcfcfb", "#1a1a19"), "panel": ("#ffffff", "#1a1a19"),
    "ink": ("#0b0b0b", "#f0efec"), "ink2": ("#52514e", "#c3c2b7"), "muted": ("#77756f", "#9a988f"),
    "grid": ("#e1e0d9", "#2c2c2a"), "axis": ("#c3c2b7", "#383835"),
    "border": ("rgba(11,11,11,0.10)", "rgba(255,255,255,0.10)"),
    "wash": ("rgba(11,11,11,0.04)", "rgba(255,255,255,0.05)"),
    "f0": ("#2a78d6", "#3987e5"), "f1": ("#eb6834", "#d95926"), "f2": ("#1baf7a", "#199e70"),
    "k-clean": ("#3d3c39", "#e6e5df"), "k-short": ("#9b9990", "#8f8d85"), "k-noframe": ("#d3d1ca", "#4a4945"),
    "r-5utr": ("#eda100", "#c98500"), "r-cds": ("#008300", "#008300"), "r-3utr": ("#e87ba4", "#d55181"),
    "r-nocds": ("#6250d6", "#9085e9"),
    "h0": ("#f0efec", "#222220"), "h1": ("#cde2fb", "#0d366b"), "h2": ("#9ec5f4", "#184f95"),
    "h3": ("#6da7ec", "#256abf"), "h4": ("#3987e5", "#3987e5"), "h5": ("#256abf", "#6da7ec"),
    "h6": ("#184f95", "#9ec5f4"), "h7": ("#0d366b", "#cde2fb"),
    "good": ("#0ca30c", "#0ca30c"), "warn": ("#fab219", "#fab219"), "crit": ("#d03b3b", "#e05252"),
    "unk": ("#898781", "#898781"),
    "lk-ink": ("#17181C", "#F4F4F1"), "lk-accent": ("#E5412D", "#FF5A43"), "lk-soft": ("#A3A8B0", "#5D626B"),
    "lk-cut": ("#FFFFFF", "#17181C"),
}
FIGURE = ("surface", "ink", "ink2", "muted", "grid", "axis", "f0", "f1", "f2", "k-clean", "k-short", "k-noframe",
          "r-5utr", "r-cds", "r-3utr", "r-nocds", "h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "warn")
C = {k: THEME[k][0] for k in FIGURE}   # the colours to draw with
FONT = "system-ui, -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"
MONO = "ui-monospace, 'SF Mono', Menlo, Consolas, 'Liberation Mono', monospace"
RC = {"font.family": "DejaVu Sans", "font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
      "ytick.labelsize": 7, "axes.edgecolor": C["axis"], "axes.linewidth": 0.6, "axes.labelcolor": C["ink2"],
      "xtick.color": C["axis"], "ytick.color": C["axis"], "xtick.labelcolor": C["muted"], "ytick.labelcolor": C["muted"],
      "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
      "text.color": C["ink"], "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": "none",
      "figure.facecolor": "none", "savefig.facecolor": "none", "savefig.edgecolor": "none", "svg.fonttype": "none",
      "grid.color": C["grid"], "grid.linewidth": 0.6,
      "lines.solid_capstyle": "round"}
HEAT_BINS = (0.01, 0.05, 0.10, 0.20, 0.30, 0.45, 0.60)   # share of the row: h0 below 0.01 ... h7 at 0.60 or more
STATE = {"pass": ("✓", "Pass"), "warn": ("!", "Warn"), "fail": ("✕", "Fail"), "unknown": ("?", "Unknown"),
         "info": ("i", "Info")}
MINUS = "−"


def esc(s):
    return html.escape(str(s), quote=True)


def tokens():
    """The CSS variables: light, then dark (the browser's choice, or data-theme on <html>)."""
    def block(i):
        return "".join(f"--{k}: {v[i]}; " for k, v in THEME.items())
    fonts = f"--sans: {FONT}; --mono: {MONO}; "
    dark = block(1) + "color-scheme: dark; "
    return (f":root {{ {block(0)}{fonts}}}\n"
            f"@media (prefers-color-scheme: dark) {{ :root:not([data-theme=\"light\"]) {{ {dark}}} }}\n"
            f":root[data-theme=\"dark\"] {{ {dark}}}\n")


def asset(name):
    """A file of this package: page.html, style.css, script.js."""
    return files(__package__).joinpath(name).read_text(encoding="utf-8")


def page(title, side, main, side_label="Flags and contents"):
    """The page: page.html with the theme's tokens, style.css and script.js in it."""
    return Template(asset("page.html")).substitute(title=esc(title), tokens=tokens(), style=asset("style.css"),
                                                   side_label=esc(side_label), side=side, main=main,
                                                   script=asset("script.js"))


def short(m):
    """A number of a path, to 0.1: 263.836364 -> 263.8, 12.0 -> 12."""
    v = f"{float(m.group()):.1f}".rstrip("0").rstrip(".")
    return "0" if v == "-0" else v


def path_data(d):
    """Path data, shorter: numbers to 0.1, no newlines, a vertical segment as V (M x y L x y2 -> Mx yVy2)."""
    d = re.sub(r"-?\d+\.\d+", short, d.replace("\n", ""))
    return re.sub(r"M (\S+) (\S+) L \1 (\S+) ?", r"M\1 \2V\3", d).strip()


def svg(fig, name, label, tips=None):
    """The figure as inline SVG: no XML header and metadata, no fixed size (the CSS scales it), colours as CSS
    variables. tips: {gid: text} becomes a <title> (a tooltip) in the element of that gid. name: unique on the page
    (the ids of the SVG come from it)."""
    import matplotlib
    buf = io.StringIO()
    with matplotlib.rc_context({"svg.hashsalt": name}):
        fig.savefig(buf, format="svg", metadata={"Date": None}, transparent=True)
    s = buf.getvalue()
    s = s[s.index("<svg"):]
    s = re.sub(r"\s*<metadata>.*?</metadata>", "", s, count=1, flags=re.S)
    s = re.sub(r'<svg ([^>]*?)width="[^"]*" height="[^"]*"', rf'<svg role="img" aria-label="{esc(label)}" \1', s,
               count=1)
    for k in FIGURE:
        s = s.replace(C[k], f"var(--{k})")
    s = s.replace("font-family: 'DejaVu Sans'", "font-family: var(--sans)")
    # path coordinates to 0.1 pt (matplotlib writes 6 decimals): a page with many bars stays small
    s = re.sub(r'\bd="([^"]*)"', lambda m: 'd="' + path_data(m.group(1)) + '"', s)
    for gid, text in (tips or {}).items():
        s = s.replace(f'<g id="{gid}">', f'<g id="{gid}"><title>{esc(text)}</title>', 1)
    # the ids of matplotlib (figure_1, axes_1, ...) repeat in each figure: prefix them, and the references to them
    s = re.sub(r'\bid="([^"]+)"', rf'id="{name}-\1"', s)
    s = re.sub(r'url\(#([^)]+)\)', rf'url(#{name}-\1)', s)
    s = re.sub(r'href="#([^"]+)"', rf'href="#{name}-\1"', s)
    left = re.findall(r"#[0-9a-fA-F]{6}\b", s)
    if left:
        raise ValueError(f"figure {name}: colours outside the theme: {sorted(set(left))}")
    return s.strip()


# ---------- numbers as text ----------

def isna(x):
    return x is None or (isinstance(x, float) and np.isnan(x)) or x is pd.NA


def pct(x, d=1):
    return "–" if isna(x) else f"{100 * x:.{d}f}%"


def f2(x):
    return "–" if isna(x) else f"{x:.2f}"


def count(x):
    return "–" if isna(x) else f"{x:,.0f}"


def nt(x):
    """A signed position: −12, 0, +3."""
    if isna(x):
        return "–"
    x = int(x)
    return f"{MINUS}{-x}" if x < 0 else f"+{x}" if x > 0 else "0"


def span(lengths):
    """[26, 27, 28, 31] -> '26–29, 31'; 'none' if empty."""
    out = [f"{r[0]}–{r[-1]}" if len(r) > 1 else f"{r[0]}" for r in qc.runs(list(lengths))]
    return ", ".join(out) or "none"


def int_list(s):
    return [int(x) for x in str(s).split(",") if x.strip()]


# ---------- page parts ----------

def icon(state):
    return f'<span class="ico {state}" aria-hidden="true">{STATE[state][0]}</span>'


def state_cell(state):
    return f'<span class="state">{icon(state)}{STATE[state][1]}</span>'


def findings(items):
    """items: (state, html) -> a list with an icon per line."""
    return '<ul class="find">' + "".join(f"<li>{icon(st)}<span>{text}</span></li>" for st, text in items) + "</ul>\n"


def tiles(items):
    """items: (label, value, note)."""
    return ('<div class="tiles">' + "".join(
        f'<div class="tile"><div class="k">{esc(k)}</div><div class="v">{esc(v)}</div><div class="s">{esc(s)}</div>'
        f"</div>" for k, v, s in items) + "</div>\n")


def html_table(df, right=(), cls="", on=()):
    """A table of df's rows as they are (strings already formatted); the columns in `right` align right, the rows at
    the positions in `on` are marked (class on)."""
    head = "".join(f'<th{" class=r" if c in right else ""}>{esc(c)}</th>' for c in df.columns)
    rows = "".join(f'<tr{" class=on" if i in on else ""}>'
                   + "".join(f'<td{" class=r" if c in right else ""}>{esc(v)}</td>' for c, v in r.items())
                   + "</tr>" for i, r in enumerate(df.to_dict("records")))
    return f'<table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table>'


def tsv_view(df, source, digits=3):
    """A disclosure with df as in its TSV (numbers rounded to `digits` decimals)."""
    shown = df.copy()
    for c in shown.columns:
        if pd.api.types.is_float_dtype(shown[c]):
            shown[c] = shown[c].map(lambda v: "NA" if isna(v) else f"{v:.{digits}f}".rstrip("0").rstrip(".")
                                    if v != int(v) else f"{int(v)}")
        else:
            shown[c] = shown[c].map(lambda v: "NA" if isna(v) else str(v))
    right = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    return (f'<details class="tbl"><summary>Table, from <code>{esc(source)}</code></summary>'
            f'<div class="scroll">{html_table(shown, right)}</div></details>\n')


def fig_block(title, chart, caption, legend="", table="", small=False):
    return (f'<figure><div class="ft">{esc(title)}</div>{legend}<div class="chart{" small" if small else ""}">'
            f"{chart}</div><figcaption>{caption}</figcaption>{table}</figure>\n")


def legend(items):
    """items: (swatch class or None, CSS background or None, label)."""
    out = []
    for cls, colour, label in items:
        style = f' style="background:{colour}"' if colour else ""
        out.append(f'<span><i class="sw{" " + cls if cls else ""}"{style}></i>{esc(label)}</span>')
    return '<div class="legend">' + "".join(out) + "</div>"


def section(sid, title, body, fold=True):
    button = '<button type="button" class="fold" aria-expanded="true">Hide</button>' if fold else ""
    return (f'<section id="{sid}" aria-labelledby="t-{sid}"><div class="sec-h"><h2 id="t-{sid}">{esc(title)}</h2>'
            f'{button}</div><div class="sec-b">\n{body}</div></section>\n')


def command(text):
    return f'<div class="cmd"><code>{esc(text)}</code><button type="button" class="copy">Copy</button></div>'

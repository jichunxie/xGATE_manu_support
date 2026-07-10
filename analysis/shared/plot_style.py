"""Shared figure style for the xGATE revision: larger axis fonts, legends placed
outside the axes, consistent palette + panel labels. Import and call apply_style()."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# study / method palettes (stable across all figures)
STUDY_PALETTE = {
    "Human CRC": "#1f77b4", "Human Prostate Cancer": "#ff7f0e",
    "Pancreas": "#2ca02c", "Senescence": "#d62728", "Liver": "#9467bd",
    "FUCCI U2OS": "#8c564b", "TS Fibroblast": "#e377c2",
}
METHOD_PALETTE = {
    "xGATE": "#d6604d", "GESECA": "#4393c3", "ORA": "#92c5de",
    "AUCell": "#999999", "scGSEA": "#762a83", "PAGODA": "#1b7837",
}
STATE_PALETTE = {"active": "#2166ac", "inactive": "#b2182b",
                 "positive": "#2166ac", "negative": "#b2182b", "partial": "#999999"}

def apply_style():
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update({
        "axes.labelsize": 16, "axes.titlesize": 15, "axes.titleweight": "bold",
        "xtick.labelsize": 13, "ytick.labelsize": 13,
        "legend.fontsize": 12, "legend.title_fontsize": 13,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
        "axes.grid": True, "grid.alpha": 0.25, "font.family": "sans-serif",
    })

def legend_outside(ax, title=None, loc="upper left", bbox=(1.02, 1.0), ncol=1, **kw):
    """Place legend to the right of the axes so it never overlaps the data."""
    return ax.legend(title=title, loc=loc, bbox_to_anchor=bbox, frameon=False,
                     borderaxespad=0.0, ncol=ncol, **kw)

def legend_below(ax, title=None, ncol=4, **kw):
    return ax.legend(title=title, loc="upper center", bbox_to_anchor=(0.5, -0.18),
                     frameon=False, ncol=ncol, **kw)

def panel_label(ax, lab, x=-0.12, y=1.08):
    # y=1.08 sits the letter clearly ABOVE a single-line centred title so the two
    # never stack. Keep titles to ONE short line (editorial wording -> caption).
    ax.text(x, y, f"({lab})", transform=ax.transAxes, fontsize=17, fontweight="bold")

def save(fig, path_noext):
    for ext in ("png", "pdf"):
        fig.savefig(f"{path_noext}.{ext}", dpi=300, bbox_inches="tight")

#!/usr/bin/env python
"""Build display items for xGATE benchmark calls to BH q<0.05.

Reads benchmark result CSVs from the repository's ``results`` and writes the
regenerated manuscript-assembly assets (``Supp_*.pdf`` / ``Figure3.jpeg`` composites and
audit CSVs). Paths are resolved through ``xgate_paths`` (config/env driven) so no absolute
paths are hard-coded; see ``configs/paths.example.yaml`` and ``docs/reproducibility.md``.
The manuscript-assembly composites are written under ``figures/manuscript_assets``
in this release (the authors' private manuscript tree is not shipped). Scientific logic and
reported values are unchanged.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd
from PIL import Image
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef
from statsmodels.stats.multitest import multipletests

from xgate_paths import ROOT, RESULTS, FIGURES, DATA


HERE = Path(__file__).resolve().parent
ROOT_DIR = Path(ROOT)
RES = Path(RESULTS)
WORK_FIG = Path(FIGURES)
# manuscript-assembly composites (Supp_*.pdf / Figure*.jpeg) — the private manuscript
# tree is not shipped; write them under the repo figures dir instead.
SUPP_FIG = Path(FIGURES) / "manuscript_assets"
MAIN_FIG = Path(FIGURES) / "manuscript_assets"
SUPP_FIG.mkdir(parents=True, exist_ok=True)

METHODS = ["xGATE", "ORA", "AUCell", "scGSEA", "GESECA", "PAGODA"]
METHOD_COLORS = {
    "xGATE": "#4C78A8",
    "ORA": "#9ecae1",
    "AUCell": "#8c8c8c",
    "scGSEA": "#9467bd",
    "GESECA": "#2ca02c",
    "PAGODA": "#e76f51",
}
DATASET_COLORS = {
    "Liver": "#4C78A8",
    "Pancreas": "#F28E2B",
    "FUCCI U2OS": "#59A14F",
    "TS Fibroblast": "#E15759",
}

# q-values are calculated over every pathway in the xGATE result file for that
# analysis, then merged into the externally labelled benchmark subset.
DATASETS = [
    ("Liver", "liver_xgate_sct.csv", "liver_sct_competing_percall.csv", "s"),
    ("Pancreas", "pancreas_xgate_bench.csv", "pancreas_sct_competing_percall.csv", "D"),
    ("FUCCI U2OS", "fucci_xgate_lognorm.csv", "fucci_lognorm_competing_percall.csv", "o"),
    ("TS Fibroblast", "ts_fibroblast_xgate_sct.csv", "ts_sct_competing_percall.csv", "^"),
]


def backup(path: Path, suffix: str) -> None:
    dest = path.with_name(path.stem + suffix + path.suffix)
    if path.exists() and not dest.exists():
        shutil.copy2(path, dest)


def truthy(values: pd.Series) -> np.ndarray:
    return values.astype(str).str.upper().isin(["TRUE", "1", "1.0"]).to_numpy()


def load_benchmark(xgate_file: str, competing_file: str):
    xg = pd.read_csv(RES / xgate_file)
    xg["pathway"] = xg["pathway"].astype(str).str.strip()
    xg["q_BH"] = multipletests(xg["p_value"], method="fdr_bh")[1]

    per = pd.read_csv(RES / competing_file)
    per.columns = [c.strip().strip('"') for c in per.columns]
    per["pathway"] = per["pathway"].astype(str).str.strip().str.strip('"')
    per["label"] = per["label"].astype(str).str.strip().str.strip('"')
    df = per.merge(xg[["pathway", "p_value", "q_BH"]], on="pathway", how="left")
    if df[["p_value", "q_BH"]].isna().any().any():
        missing = df.loc[df.q_BH.isna(), "pathway"].tolist()
        raise RuntimeError(f"Missing xGATE results for {missing}")

    truth = df.label.eq("positive").to_numpy()
    calls = {"xGATE": df.q_BH.lt(0.05).to_numpy()}
    scores = {"xGATE": (1.0 - df.p_value).to_numpy()}
    for method, active_col in [
        ("ORA", "ORA_active"),
        ("AUCell", "AUCell_active"),
        ("scGSEA", "scGSEA_active"),
        ("GESECA", "GESECA_active"),
        ("PAGODA", "PAGODA_active"),
    ]:
        calls[method] = truthy(df[active_col])
        score_col = active_col.replace("_active", "_score")
        scores[method] = pd.to_numeric(df[score_col], errors="coerce").to_numpy()
    return df, truth, calls, scores


def metrics(y: np.ndarray, call: np.ndarray, score: np.ndarray) -> dict:
    y = np.asarray(y, bool)
    call = np.asarray(call, bool)
    tp = int((y & call).sum())
    fp = int((~y & call).sum())
    fn = int((y & ~call).sum())
    tn = int((~y & ~call).sum())
    precision = tp / (tp + fp) if tp + fp else np.nan
    recall = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    specificity = tn / (tn + fp) if tn + fp else np.nan
    mcc = matthews_corrcoef(y, call) if len(set(y)) > 1 and len(set(call)) > 1 else np.nan
    score = np.asarray(score, float)
    aucpr = average_precision_score(y, score) if np.isfinite(score).all() else np.nan
    return {
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "specificity": specificity,
        "accuracy": (tp + tn) / len(y),
        "MCC": mcc,
        "AUCPR": aucpr,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
    }


def render_call_table(name: str, df: pd.DataFrame, calls: dict, output: Path) -> None:
    positives = df.loc[df.label.eq("positive"), "pathway"].tolist()
    negatives = df.loc[df.label.eq("negative"), "pathway"].tolist()
    callmap = {m: dict(zip(df.pathway, calls[m])) for m in METHODS}
    max_len = max(len(positives), len(negatives))
    fig, axes = plt.subplots(1, 2, figsize=(13, 0.47 * max_len + 3.2))
    for ax, title, pathways, expected in [
        (axes[0], "Active / positive", positives, True),
        (axes[1], "Inactive / negative", negatives, False),
    ]:
        ax.set_xlim(0, len(METHODS))
        ax.set_ylim(0, max_len)
        ax.invert_yaxis()
        ax.axis("off")
        ax.text(
            len(METHODS) / 2,
            -2.6,
            title,
            ha="center",
            va="center",
            fontsize=15,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", fc="#cfe2f3", ec="#9fc5e8"),
            clip_on=False,
        )
        for j, method in enumerate(METHODS):
            ax.text(j + 0.5, -0.1, method, ha="left", va="bottom", fontsize=12,
                    fontweight="bold", rotation=45)
        for i, pathway in enumerate(pathways):
            ax.text(-0.1, i + 0.5, pathway[:42], ha="right", va="center", fontsize=11)
            for j, method in enumerate(METHODS):
                correct = bool(callmap[method][pathway]) == expected
                ax.text(j + 0.5, i + 0.5, "✓" if correct else "✗", ha="center", va="center",
                        fontsize=17, color="#2e7d32" if correct else "#c62828")
    fig.suptitle(f"{name} benchmark: correct (✓) / incorrect (✗) by method",
                 fontsize=15, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_benchmark_assets() -> pd.DataFrame:
    rows = []
    loaded = {}
    output_names = {
        "Liver": "Supp_liver_benchmark_pathway_activity.pdf",
        "Pancreas": "Supp_pancreas_benchmark_pathway_activity.pdf",
        "FUCCI U2OS": "Supp_fucci_benchmark.pdf",
        "TS Fibroblast": "Supp_ts_benchmark.pdf",
    }
    for name, xgate_file, competing_file, marker in DATASETS:
        df, truth, calls, scores = load_benchmark(xgate_file, competing_file)
        loaded[name] = (df, truth, calls, scores, marker)
        for method in METHODS:
            row = metrics(truth, calls[method], scores[method])
            row.update(dataset=name, method=method)
            rows.append(row)
        out = SUPP_FIG / output_names[name]
        backup(out, ".rawp_backup")
        render_call_table(name, df, calls, out)

    metric_df = pd.DataFrame(rows)[
        ["dataset", "method", "precision", "recall", "F1", "specificity", "accuracy",
         "MCC", "AUCPR", "TP", "FP", "FN", "TN"]
    ]
    metric_df.to_csv(HERE / "fig3_benchmark_metrics_bh.csv", index=False)
    build_figure3(metric_df, loaded)
    return metric_df


def panel_label(ax, label: str, x=-0.12, y=1.05):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top")


def build_figure3(metric_df: pd.DataFrame, loaded: dict) -> None:
    old_path = MAIN_FIG / "Figure3.jpeg"
    backup(old_path, ".rawp_backup")
    old = Image.open(old_path).convert("RGB")
    split_y = 1180
    bottom = old.crop((0, split_y, old.width, old.height))

    fig = plt.figure(figsize=(33.25, 3.93), dpi=300)
    outer = GridSpec(1, 2, figure=fig, width_ratios=[0.37, 0.63], wspace=0.16)
    left = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[0.47, 0.53], wspace=0.28)
    right = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1], wspace=0.28)
    ax_pr = fig.add_subplot(left[0])
    ax_mean = fig.add_subplot(left[1])

    markers = {name: item[4] for name, item in loaded.items()}
    for _, row in metric_df.iterrows():
        ax_pr.scatter(row.recall, row.precision, s=105,
                      color=METHOD_COLORS[row.method], marker=markers[row.dataset],
                      edgecolor="0.25", linewidth=0.5, alpha=0.9, zorder=3)
    rr = np.linspace(0.01, 1, 200)
    for f1 in (0.5, 0.7, 0.9):
        pp = f1 * rr / (2 * rr - f1)
        pp[(pp < 0) | (pp > 1)] = np.nan
        ax_pr.plot(rr, pp, ":", color="0.75", lw=0.8)
        ax_pr.text(1.01, f1 / (2 - f1), f"F1={f1:.1f}", fontsize=8, va="center")
    ax_pr.set_xlim(0, 1.14)
    ax_pr.set_ylim(0, 1.05)
    ax_pr.set_xlabel("Recall")
    ax_pr.set_ylabel("Precision")
    ax_pr.set_title("Precision and recall", fontweight="bold")
    panel_label(ax_pr, "a", x=-0.16, y=1.12)

    method_handles = [mlines.Line2D([], [], color=METHOD_COLORS[m], marker="o", ls="", ms=5, label=m)
                      for m in METHODS]
    dataset_handles = [mlines.Line2D([], [], color="0.25", marker=markers[d], ls="", ms=5, label=d)
                       for d in markers]
    leg_method = fig.legend(handles=method_handles, title="Method", frameon=False, fontsize=5.5,
                            title_fontsize=6, loc="center", bbox_to_anchor=(0.188, 0.68))
    fig.legend(handles=dataset_handles, title="Dataset marker", frameon=False, fontsize=5.5,
               title_fontsize=6, loc="center", bbox_to_anchor=(0.188, 0.33))
    fig.add_artist(leg_method)

    mean = metric_df.groupby("method", sort=False)[
        ["precision", "recall", "F1", "specificity", "MCC", "AUCPR"]
    ].mean().reindex(METHODS)
    x = np.arange(len(mean.columns))
    width = 0.12
    for i, method in enumerate(METHODS):
        ax_mean.bar(x + (i - 2.5) * width, mean.loc[method], width,
                    color=METHOD_COLORS[method], label=method)
    ax_mean.axhline(0, color="0.4", lw=0.6)
    ax_mean.set_xticks(x, mean.columns, rotation=20, ha="right")
    ax_mean.set_ylim(-0.25, 1.05)
    ax_mean.set_ylabel("score (mean of 4 datasets)")
    ax_mean.set_title("Mean benchmark metrics", fontweight="bold")

    for col, ax in zip(["specificity", "MCC", "AUCPR"], [fig.add_subplot(right[i]) for i in range(3)]):
        pivot = metric_df.pivot(index="method", columns="dataset", values=col).reindex(METHODS)
        xx = np.arange(len(METHODS))
        w = 0.18
        for j, dataset in enumerate(DATASET_COLORS):
            ax.bar(xx + (j - 1.5) * w, pivot[dataset], w, color=DATASET_COLORS[dataset], label=dataset)
        ax.axhline(0, color="0.4", lw=0.6)
        ax.set_xticks(xx, METHODS, rotation=25, ha="right")
        ax.set_title(col, fontweight="bold")
        ax.set_ylim((-0.6, 1.05) if col == "MCC" else (0, 1.05))
    right_axes = fig.axes[-3:]
    panel_label(right_axes[0], "b", x=-0.16, y=1.12)
    handles = [mlines.Line2D([], [], color=c, marker="s", ls="", ms=8, label=d)
               for d, c in DATASET_COLORS.items()]
    right_axes[-1].legend(handles=handles, title="Dataset", frameon=False,
                          loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8)

    for ax in fig.axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.18)
    fig.subplots_adjust(left=0.035, right=0.97, top=0.88, bottom=0.22)

    tmp = HERE / "Figure3_bh_top.png"
    fig.savefig(tmp, dpi=300, facecolor="white")
    plt.close(fig)
    top = Image.open(tmp).convert("RGB").resize((old.width, split_y), Image.Resampling.LANCZOS)
    combined = Image.new("RGB", old.size, "white")
    combined.paste(top, (0, 0))
    combined.paste(bottom, (0, split_y))
    combined.save(old_path, quality=95, dpi=(300, 300))


def build_real_bh_by_figure() -> pd.DataFrame:
    rows = []
    pvalue_files = [
        ("Pancreas", "supp_figS11_pancreas_ctrl_focused.csv"),
        ("Liver", "liver_xgate_sct.csv"),
        ("FUCCI U2OS", "fucci_xgate_lognorm.csv"),
        ("TS Fibroblast", "ts_fibroblast_xgate_sct.csv"),
    ]
    for name, filename in pvalue_files:
        df = pd.read_csv(RES / filename).dropna(subset=["p_value"]).drop_duplicates("pathway")
        q_bh = multipletests(df.p_value, method="fdr_bh")[1]
        q_by = multipletests(df.p_value, method="fdr_by")[1]
        rows.append({
            "dataset": name,
            "n_tested": len(df),
            "n_sig_BH": int((q_bh < 0.05).sum()),
            "n_sig_BY": int((q_by < 0.05).sum()),
        })
    summary = pd.DataFrame(rows)
    summary["fraction_sig_BH"] = summary.n_sig_BH / summary.n_tested
    summary["fraction_sig_BY"] = summary.n_sig_BY / summary.n_tested
    summary.to_csv(HERE / "bh_by_real_summary.csv", index=False)

    source = Image.open(WORK_FIG / "fig3c_jaccard_pvalue.png").convert("RGB")
    top = source.crop((0, 0, source.width, source.height // 2 + 30))
    fig = plt.figure(figsize=(20, 10))
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1.05, 1.0], hspace=0.24, wspace=0.28)
    ax_top = fig.add_subplot(gs[0, :])
    ax_top.imshow(top)
    ax_top.axis("off")
    axes = [fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2:])]
    for ax, method, count_col, frac_col, label in [
        (axes[0], "BH", "n_sig_BH", "fraction_sig_BH", "e"),
        (axes[1], "BY", "n_sig_BY", "fraction_sig_BY", "f"),
    ]:
        x = np.arange(len(summary))
        bars = ax.bar(x, summary[count_col], color=[DATASET_COLORS[d] for d in summary.dataset],
                      edgecolor="0.25", linewidth=0.7)
        ax.set_xticks(x, summary.dataset, rotation=15, ha="right")
        ax.set_ylabel("pathways significant at q < 0.05")
        ax.set_title(f"{method}: significant pathways in the real datasets", fontweight="bold")
        ax.set_ylim(0, max(summary.n_sig_BH.max(), 1) * 1.23)
        for bar, n_sig, n_tested, frac in zip(
            bars, summary[count_col], summary.n_tested, summary[frac_col]
        ):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.7,
                    f"{n_sig}/{n_tested} ({100 * frac:.0f}%)", ha="center", va="bottom", fontsize=10)
        panel_label(ax, label, x=-0.05, y=1.06)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.2)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.98, bottom=0.09)
    out_pdf = SUPP_FIG / "Supp_jaccard_pvalue.pdf"
    backup(out_pdf, ".simulation_backup")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(HERE / "Supp_jaccard_pvalue_real_bh_by.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    return summary


def build_complexity_figure() -> tuple[float, float, float]:
    base = pd.read_csv(ROOT_DIR / "computational_benchmark_summary.csv")
    additions = [pd.read_csv(RES / "computational_benchmark_summary_additions.csv")]
    additions.append(pd.read_csv(Path(DATA) / "ts_fibroblast/adj_matrix_ts_fibroblast_pooled_timing.csv"))
    raw = pd.concat([base] + additions, ignore_index=True)
    raw["dataset_label"] = raw.dataset_label.fillna(raw.dataset)
    vertices = raw.n_vertices.fillna(raw.n_genes).astype(float)
    edges = raw.n_edges.astype(float)
    df = pd.DataFrame({
        "dataset": raw.dataset_label,
        "candidate_pairs": vertices * (vertices - 1) / 2,
        "edges": edges,
        "density": edges / (vertices * (vertices - 1) / 2),
        "cpu": pd.to_numeric(raw.cpu_total_s, errors="coerce"),
        "memory": pd.to_numeric(raw.peak_rss_mb, errors="coerce") / 1024,
    }).dropna(subset=["candidate_pairs", "edges", "cpu"])
    rho_cpu = spearmanr(df.candidate_pairs, df.cpu).statistic
    rho_mem = spearmanr(df.candidate_pairs, df.memory, nan_policy="omit").statistic
    rho_density = spearmanr(df.density, df.cpu).statistic

    order = ["Human CRC", "Human Prostate Cancer", "Pancreas", "Senescence", "Liver",
             "FUCCI U2OS", "TS Fibroblast"]
    palette = dict(zip(order, sns.color_palette("tab10", n_colors=len(order))))
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))

    def scatter(ax, x, y, xlabel, ylabel, title):
        sns.scatterplot(data=df, x=x, y=y, hue="dataset", style="dataset", hue_order=order,
                        style_order=order, palette=palette, s=95, alpha=0.85, ax=ax, legend=False)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontweight="bold", pad=13)

    scatter(axes[0, 0], "candidate_pairs", "cpu", r"Candidate gene pairs $|V|(|V|-1)/2$",
            "CPU time (s)", rf"CPU time vs candidate pairs; Spearman $\rho={rho_cpu:.2f}$")
    slope, intercept = np.polyfit(np.log10(df.candidate_pairs), np.log10(df.cpu), 1)
    xx = np.array([df.candidate_pairs.min(), df.candidate_pairs.max()])
    axes[0, 0].plot(xx, 10 ** intercept * xx ** slope, "k--", lw=1.2, zorder=0)
    scatter(axes[0, 1], "candidate_pairs", "memory", r"Candidate gene pairs $|V|(|V|-1)/2$",
            "Peak memory (GB)", rf"Peak memory vs candidate pairs; Spearman $\rho={rho_mem:.2f}$")
    scatter(axes[1, 0], "density", "cpu", "Realized edge density", "CPU time (s)",
            rf"CPU time vs edge density; Spearman $\rho={rho_density:.2f}$ (n.s.)")
    sns.boxplot(data=df, x="dataset", y="density", order=order, hue="dataset", palette=palette,
                legend=False, showfliers=False, ax=axes[1, 1])
    sns.stripplot(data=df, x="dataset", y="density", order=order, color="0.25", size=3,
                  alpha=0.5, ax=axes[1, 1])
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("Study")
    axes[1, 1].set_ylabel(r"Edge density $|E|/[|V|(|V|-1)/2]$")
    axes[1, 1].set_title("Edge density by study", fontweight="bold", pad=13)
    plt.setp(axes[1, 1].get_xticklabels(), rotation=30, ha="right")
    for ax, label in zip(axes.ravel(), "abcd"):
        panel_label(ax, label, x=-0.13, y=1.05)
        ax.spines[["top", "right"]].set_visible(False)
    handles = [mlines.Line2D([], [], color=palette[d], marker="o", ls="", ms=8, label=d)
               for d in order]
    fig.legend(handles=handles, title="Study", loc="center left", bbox_to_anchor=(0.91, 0.5),
               frameon=False)
    fig.tight_layout(rect=[0, 0, 0.9, 1])
    out = SUPP_FIG / "Supp_complexity_analysis.pdf"
    backup(out, ".slope_backup")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(HERE / "Supp_complexity_analysis_spearman.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    return rho_cpu, rho_mem, rho_density


def main() -> None:
    metric_df = build_benchmark_assets()
    bh_by = build_real_bh_by_figure()
    rhos = build_complexity_figure()
    xgate = metric_df.loc[metric_df.method.eq("xGATE"), ["dataset", "F1", "MCC", "AUCPR"]]
    print("\nBH-controlled xGATE benchmark metrics")
    print(xgate.round(3).to_string(index=False))
    print(f"mean F1 = {xgate.F1.mean():.3f}")
    print("\nReal-data BH/BY summary")
    print(bh_by.to_string(index=False))
    print(f"\nComplexity Spearman rho: CPU={rhos[0]:.3f}, memory={rhos[1]:.3f}, density={rhos[2]:.3f}")


if __name__ == "__main__":
    main()

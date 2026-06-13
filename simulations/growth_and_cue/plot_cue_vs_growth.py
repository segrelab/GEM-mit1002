from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.rcParams.update(
    {
        "font.size": 11,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]
IN_PATH = FILE_PATH / "results"
OUT_PATH = FILE_PATH / "figures"
OUT_PATH.mkdir(exist_ok=True)

# Import plot_styles.py from the root of the repo
sys.path.append(str(REPO_ROOT))
from plot_styles import summer_colors

# SMF-figure palette, extended with two harmonising tones for the 7 categories.
ENTRY_POINT_COLORS = {
    "Glycolysis": summer_colors["light_blue"],
    "Pyruvate": summer_colors["teal"],
    "Acetyl-CoA": summer_colors["yellow"],
    "TCA — α-KG": summer_colors["pink"],
    "TCA — OAA": summer_colors["dark_pink"],
    "TCA — Succinyl-CoA": summer_colors["green"],
    "Aromatic catabolism": "#9B6A9E",  # muted mauve (added)
}
ENTRY_POINT_ORDER = list(ENTRY_POINT_COLORS)


def main():
    # Load the Growth and CUE results
    summary_df = pd.read_csv(IN_PATH / "growth_and_cue.csv", index_col=0)
    # Plot the bar charts of growth rate and CUE, coloured by entry point into central metabolism
    plot_growth_cue(summary_df, OUT_PATH / "growth_and_cue.png")
    # Plot the correlation between growth rate and CUE across substrates
    plot_growth_cue_correlation(summary_df, OUT_PATH / "growth_vs_cue_scatter.png")


def plot_growth_cue(summary_df: pd.DataFrame, out_path: Path) -> None:
    df = summary_df.sort_values("growth_rate", ascending=False).reset_index()
    colors = [ENTRY_POINT_COLORS.get(c, "#9E9E9E") for c in df["entry_point"]]
    x = np.arange(len(df))

    fig, (ax_g, ax_c) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

    ax_g.bar(x, df["growth_rate"], color=colors, edgecolor="white", linewidth=0.6)
    ax_g.set_ylabel("Growth rate (h⁻¹)")
    ax_g.set_title(
        "MIT1002 growth rate and carbon-use efficiency across substrates",
        fontsize=12,
        pad=8,
    )
    ax_g.margins(x=0.01)

    ax_c.bar(x, df["cue"], color=colors, edgecolor="white", linewidth=0.6)
    ax_c.set_ylabel("Carbon-use efficiency")
    ax_c.set_ylim(0, 1)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(df["substrate"], rotation=45, ha="right", fontsize=9)

    handles = [
        mpatches.Patch(facecolor=ENTRY_POINT_COLORS[lbl], edgecolor="white", label=lbl)
        for lbl in ENTRY_POINT_ORDER
        if lbl in df["entry_point"].values
    ]
    ax_g.legend(
        handles=handles,
        title="Entry point into\ncentral metabolism",
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
        title_fontsize=9.5,
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_growth_cue_correlation(summary_df: pd.DataFrame, out_path: Path) -> None:
    df = summary_df.reset_index()
    colors = [ENTRY_POINT_COLORS.get(c, "#9E9E9E") for c in df["entry_point"]]
    x = df["growth_rate"].values.astype(float)
    y = df["cue"].values.astype(float)

    r, p = stats.pearsonr(x, y)
    rho, p_s = stats.spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.scatter(x, y, c=colors, s=120, edgecolor="white", linewidth=0.6, zorder=3)

    # Least-squares fit line
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, slope * xs + intercept, color="#555555", ls="--", lw=1.2, zorder=2)

    # Substrate labels. Default is offset to the right of the point; a few
    # near-coincident points get a manual override so labels don't overlap.
    label_offsets = {  # substrate: (dx, dy, ha, va)
        "Galactose": (0, 9, "center", "bottom"),  # sits almost on Glucose
    }
    for _, row in df.iterrows():
        dx, dy, ha, va = label_offsets.get(row["substrate"], (5, 0, "left", "center"))
        ax.annotate(
            row["substrate"],
            (row["growth_rate"], row["cue"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7,
            ha=ha,
            va=va,
            color="#555555",
        )

    ax.set_xlabel("Growth rate (h⁻¹)")
    ax.set_ylabel("Carbon-use efficiency")
    ax.set_title(
        "Growth rate vs. carbon-use efficiency across substrates", fontsize=12, pad=8
    )

    # Correlation metrics
    txt = (
        f"Pearson  r = {r:.2f}  (p = {p:.2g})\n"
        f"Spearman ρ = {rho:.2f}  (p = {p_s:.2g})"
    )
    ax.text(
        0.025,
        0.025,
        txt,
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="left",
        linespacing=1.4,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc", alpha=0.9),
    )

    # Legend: entry-point colors + the fit line
    handles = [
        mpatches.Patch(facecolor=ENTRY_POINT_COLORS[lbl], edgecolor="white", label=lbl)
        for lbl in ENTRY_POINT_ORDER
        if lbl in df["entry_point"].values
    ]
    handles.append(
        Line2D([0], [0], color="#555555", ls="--", lw=1.2, label="linear fit")
    )
    ax.legend(
        handles=handles,
        title="Entry point into\ncentral metabolism",
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
        title_fontsize=9.5,
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}  (Pearson r={r:.2f}, p={p:.2g})")


if __name__ == "__main__":
    main()

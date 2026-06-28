"""Plot the cost of forcing flux through the ED pathway.

Reproduces the growth-rate, CUE, BGE, and acetate-secretion figures from
plot_growth_vs_ed_use.ipynb, styled with the repo's plot_styles module and the
"summer" color palette.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Import plot_styles.py from the root of the repo
sys.path.append(str(REPO_ROOT))
from plot_styles import set_plot_style, summer_colors  # noqa: E402

# Color each O2 level with a summer-palette color
O2_COLORS = {20: summer_colors["teal"], 1000: summer_colors["pink"]}


def plot_metric(df, metric, ylabel, title, filename, ylim=None):
    fig, ax = plt.subplots()
    for o2_level in sorted(df["o2"].unique()):
        subset = df[df["o2"] == o2_level].sort_values("ed_percent")
        ax.plot(
            subset["ed_percent"],
            subset[metric],
            marker="o",
            color=O2_COLORS.get(o2_level, "gray"),
            label=f"O$_2$: {o2_level}",
        )
    ax.set_xlabel("Percent of Glycolytic Flux through ED")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend()
    set_plot_style(ax)
    fig.tight_layout()
    fig.savefig(OUT_PATH / filename, dpi=150)
    plt.close(fig)


def main():
    df = pd.read_csv(OUT_PATH / "growth_vs_ed_use_results.csv")

    # Growth rate (the headline tradeoff): default autoscaled + a 0-anchored view
    plot_metric(
        df, "growth_rate", "Growth Rate (h$^{-1}$)", "Growth Rate vs ED Use",
        "growth_rate_vs_ed_use.png",
    )
    plot_metric(
        df, "growth_rate", "Growth Rate (h$^{-1}$)", "Growth Rate vs ED Use",
        "growth_rate_vs_ed_use_full_range.png", ylim=(0, 0.85),
    )
    plot_metric(
        df, "cue", "Carbon Use Efficiency", "CUE vs ED Use",
        "cue_vs_ed_use.png", ylim=(0, 1),
    )
    plot_metric(
        df, "bge", "Bacterial Growth Efficiency (G / G+R)", "BGE vs ED Use",
        "bge_vs_ed_use.png", ylim=(0, 1),
    )
    plot_metric(
        df, "acetate_secretion", "Acetate Secretion (mmol/gDW/hr)",
        "Acetate Secretion vs ED Use", "acetate_secretion_vs_ed_use.png",
    )
    print("Saved figures to", OUT_PATH)


if __name__ == "__main__":
    main()

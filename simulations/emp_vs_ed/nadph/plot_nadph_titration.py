"""Plot the NADPH titration results.

Figure 1 (the money shot): percent of glycolytic flux through ED as the forced
NADPH demand rises, one line per O2 level, faceted by substrate. On glucose ED
should climb from ~0 with redox demand; on galacturonate ED is pinned high
because it is the only route in.

Figure 2: how the NADPH demand is actually met on glucose -- ED versus the
competing valves (oxPPP, NADP-IDH, malic enzyme, transhydrogenase) -- faceted
by O2 level.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Import the shared plot styles from tools/
sys.path.append(str(REPO_ROOT))
from tools.plot_styles import set_plot_style, summer_colors  # noqa: E402

# NADPH source columns in the order/colors we want to show them
SOURCE_COLS = [
    "ed_flux",
    "6PGDH / oxPPP (rxn01115)",
    "NADP-IDH (rxn01387)",
    "Malic enzyme (rxn00161)",
    "Transhydrogenase (rxn09295)",
]
SOURCE_LABELS = {
    "ed_flux": "ED pathway",
    "6PGDH / oxPPP (rxn01115)": "oxidative PPP",
    "NADP-IDH (rxn01387)": "NADP-IDH",
    "Malic enzyme (rxn00161)": "malic enzyme",
    "Transhydrogenase (rxn09295)": "transhydrogenase",
}
SOURCE_COLORS = [
    summer_colors["pink"],
    summer_colors["teal"],
    summer_colors["yellow"],
    summer_colors["green"],
    summer_colors["dark_tan"],
]

# Colors for the O2 levels in figure 1
O2_COLORS = {20: summer_colors["teal"], 1000: summer_colors["pink"]}


def load_results():
    df = pd.read_csv(OUT_PATH / "nadph_titration_results.csv")
    # Keep only feasible solves
    return df[df["status"] == "optimal"].copy()


def plot_percent_ed(df):
    """Figure 1: percent ED flux vs forced NADPH demand, faceted by substrate."""
    substrates = sorted(df["substrate"].unique())
    fig, axes = plt.subplots(
        1, len(substrates), figsize=(5 * len(substrates), 5), sharey=True
    )
    if len(substrates) == 1:
        axes = [axes]

    for ax, substrate in zip(axes, substrates):
        sub = df[df["substrate"] == substrate]
        for o2_level, grp in sub.groupby("o2_level"):
            grp = grp.sort_values("nadph_drain")
            ax.plot(
                grp["nadph_drain"],
                grp["percent_ed_flux"] * 100,
                marker="o",
                markersize=3,
                color=O2_COLORS.get(o2_level, "gray"),
                label=f"O$_2$ = {o2_level}",
            )
        ax.set_title(substrate, color="gray")
        ax.set_xlabel("Forced NADPH demand (mmol gDW$^{-1}$ h$^{-1}$)")
        ax.set_ylim(-2, 102)
        ax.legend()
        set_plot_style(ax)

    axes[0].set_ylabel("Percent of glycolytic flux through ED")
    fig.tight_layout()
    fig.savefig(OUT_PATH / "percent_ed_vs_nadph_demand.png", dpi=150)
    plt.close(fig)


def plot_source_breakdown(df, substrate="Glucose"):
    """Figure 2: which NADPH valves carry the demand, faceted by O2 level."""
    sub = df[df["substrate"] == substrate]
    o2_levels = sorted(sub["o2_level"].unique())
    fig, axes = plt.subplots(
        1, len(o2_levels), figsize=(5 * len(o2_levels), 5), sharey=True
    )
    if len(o2_levels) == 1:
        axes = [axes]

    for ax, o2_level in zip(axes, o2_levels):
        grp = sub[sub["o2_level"] == o2_level].sort_values("nadph_drain")
        for col, color in zip(SOURCE_COLS, SOURCE_COLORS):
            if col not in grp.columns:
                continue
            ax.plot(
                grp["nadph_drain"],
                grp[col],
                marker="o",
                markersize=3,
                color=color,
                label=SOURCE_LABELS.get(col, col),
            )
        ax.set_title(f"{substrate}, O$_2$ = {o2_level}", color="gray")
        ax.set_xlabel("Forced NADPH demand (mmol gDW$^{-1}$ h$^{-1}$)")
        ax.legend()
        set_plot_style(ax)

    axes[0].set_ylabel("Reaction flux (mmol gDW$^{-1}$ h$^{-1}$)")
    fig.tight_layout()
    fig.savefig(OUT_PATH / "nadph_source_breakdown.png", dpi=150)
    plt.close(fig)


def main():
    df = load_results()
    plot_percent_ed(df)
    plot_source_breakdown(df, substrate="Glucose")
    print("Saved figures to", OUT_PATH)


if __name__ == "__main__":
    main()

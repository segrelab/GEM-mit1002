"""Plot the oxygen cost of growing on phenol.

Figure 1: growth rate vs O2 supply for phenol vs glucose vs acetate.
Figure 2: O2 budget per substrate carbon, split into catabolic (oxygenase) and
respiratory O2.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

sys.path.append(str(REPO_ROOT))
from plot_styles import set_plot_style, summer_colors  # noqa: E402

SUBSTRATE_COLORS = {
    "Phenol": summer_colors["pink"],
    "Glucose": summer_colors["teal"],
    "Acetate": summer_colors["yellow"],
}


def plot_growth_vs_o2(df):
    fig, ax = plt.subplots()
    for substrate in ["Glucose", "Acetate", "Phenol"]:
        sub = df[df["substrate"] == substrate].sort_values("o2_supply")
        ax.plot(
            sub["o2_supply"], sub["growth_rate"], marker="o", markersize=3,
            color=SUBSTRATE_COLORS.get(substrate, "gray"), label=substrate,
        )
    ax.set_xlabel("O$_2$ supply (mmol gDW$^{-1}$ h$^{-1}$)")
    ax.set_ylabel("Growth rate (h$^{-1}$)")
    ax.set_title("Growth vs O$_2$ supply (fixed carbon uptake)")
    ax.legend()
    set_plot_style(ax)
    fig.tight_layout()
    fig.savefig(OUT_PATH / "growth_vs_o2.png", dpi=150)
    plt.close(fig)


def plot_o2_budget(df):
    df = df.set_index("substrate").loc[["Glucose", "Acetate", "Phenol"]]
    x = np.arange(len(df))
    fig, ax = plt.subplots()
    resp = df["respiratory_o2_per_c"].values
    oxy = df["oxygenase_o2_per_c"].values
    ax.bar(x, resp, color=summer_colors["teal"], label="Respiratory (terminal oxidase)")
    ax.bar(x, oxy, bottom=resp, color=summer_colors["pink"],
           label="Oxygenase (catabolic co-substrate)")
    ax.set_xticks(x)
    ax.set_xticklabels(df.index)
    ax.set_xlabel("Substrate")
    ax.set_ylabel("O$_2$ consumed per substrate C (mol/mol)")
    ax.set_title("O$_2$ budget per carbon")
    ax.legend()
    set_plot_style(ax)
    fig.tight_layout()
    fig.savefig(OUT_PATH / "o2_budget_per_carbon.png", dpi=150)
    plt.close(fig)


def main():
    plot_growth_vs_o2(pd.read_csv(OUT_PATH / "growth_vs_o2.csv"))
    plot_o2_budget(pd.read_csv(OUT_PATH / "o2_budget.csv"))
    print("Saved figures to", OUT_PATH)


if __name__ == "__main__":
    main()

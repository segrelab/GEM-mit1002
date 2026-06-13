from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

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

# Color palette for exchange metabolites
# The "Summer" color palette with a few extra colors to avoid repeats
# Remove colors that I hard set below
EX_PALETTE = [
    v
    for k, v in summer_colors.items()
    if k not in {"yellow", "light_blue", "dark_tan", "light_tan"}
]
EX_PALETTE += [
    "#9B6A9E",  # mauve
    "#3F6F8C",  # steel blue
    "#C98B8B",  # dusty rose
    "#8C4F6B",  # plum
    "#6E7B8B",  # slate
    "#BFB13A",  # brass
]
# Define an eye-catching color for ammonium, since that is the primary result
NH3_COLOR = summer_colors["yellow"]
# Define muted colors for metabolites that are unimportant
H_COLOR = summer_colors["light_blue"]
CO2_COLOR = summer_colors["dark_tan"]
H2O_COLOR = summer_colors["light_tan"]
# Define a single color for the carbon source
# So there aren't unique colros for each bar
# OK to be the same color as CO2, since the carbon source is ony ever taken up
# And CO2 is only ever released
CARBON_SOURCE_COLOR = summer_colors["dark_tan"]  # all substrate carbon sources
OTHER_COLOR = "#bdbdbd"  # grey — collapsed trace metabolites


def main():
    # Load the exchange fluxes
    ex_df = pd.read_csv(IN_PATH / "exchange_fluxes.csv", index_col=0)

    # Load the substrate panel to get the names of the carbon sources
    substrate_df = pd.read_csv(IN_PATH / "substrate_panel.csv")
    carbon_source_names = substrate_df["name_in_model"].tolist()

    # Plot
    plot_exchange_stacks(ex_df, carbon_source_names, OUT_PATH)


def plot_exchange_stacks(ex_df, carbon_source_names, out_dir):
    """Two stacked-bar charts, both drawn upward, with shared metabolite colours.

    Uptake fluxes are negated so they read as positive bars. All substrate
    carbon sources are merged into one 'Carbon source' segment in the uptake
    chart (each appears on a single bar). Byproducts/co-substrates keep
    consistent colours across both charts.
    """
    # Separate uptake and exudation
    uptake = -ex_df.clip(upper=0)  # negative fluxes -> positive uptake
    exud = ex_df.clip(lower=0)  # positive fluxes -> exudation

    # Drop zero-only columns (e.g. if no exudation of a metabolite occurs)
    uptake = uptake.loc[:, (uptake != 0).any(axis=0)]
    exud = exud.loc[:, (exud != 0).any(axis=0)]

    # Merge all substrate carbon sources into one column in the uptake chart
    # So they don't have unique colours and the chart is less busy
    uptake = _merge_carbon_sources(uptake, carbon_source_names)

    # Drop H2O from the exudation chart
    exud = exud.drop(columns=["H2O"], errors="ignore")

    # Shared colour map over the "regular" metabolites (everything except the
    # fixed-colour Carbon source / Other), ordered by total magnitude so shared
    # metabolites get the same colour in both charts.
    special = {"H2O", "H+", "CO2", "NH3", "Carbon source", "Other"}
    regular = (set(uptake.columns) | set(exud.columns)) - special

    def total_mag(c):
        m = 0.0
        if c in uptake.columns:
            m += uptake[c].abs().sum()
        if c in exud.columns:
            m += exud[c].abs().sum()
        return m

    ordered_reg = sorted(regular, key=total_mag, reverse=True)
    colors = {c: EX_PALETTE[i % len(EX_PALETTE)] for i, c in enumerate(ordered_reg)}
    # Define the special colors
    colors["H+"] = H_COLOR
    colors["H2O"] = H2O_COLOR
    colors["CO2"] = CO2_COLOR
    colors["NH3"] = NH3_COLOR
    colors["Carbon source"] = CARBON_SOURCE_COLOR
    colors["Other"] = OTHER_COLOR

    up_cols = _ordered_cols(uptake, pin_first="Carbon source", pin_last="Other")
    ex_cols = _ordered_cols(exud, pin_last="Other")
    _stacked_bar(
        uptake[up_cols],
        colors,
        "MIT1002 uptake fluxes across substrates",
        "Uptake flux (mmol gDW⁻¹ h⁻¹)",
        out_dir / "uptake_fluxes.png",
    )
    _stacked_bar(
        exud[ex_cols],
        colors,
        "MIT1002 exudation fluxes across substrates",
        "Exudation flux (mmol gDW⁻¹ h⁻¹)",
        out_dir / "exudation_fluxes.png",
    )


def _merge_carbon_sources(df, carbon_source_names):
    """Collapse all substrate-carbon-source columns into one 'Carbon source'."""
    cs = [c for c in df.columns if c in carbon_source_names]
    rest = [c for c in df.columns if c not in carbon_source_names]
    out = df[rest].copy()
    if cs:
        out["Carbon source"] = df[cs].sum(axis=1)
    return out


def _ordered_cols(df, pin_first=None, pin_last=None):
    """Column order: pinned-first, then regular by descending magnitude, pinned-last."""
    skip = {pin_first, pin_last}
    middle = sorted(
        (c for c in df.columns if c not in skip),
        key=lambda c: df[c].abs().sum(),
        reverse=True,
    )
    head = [pin_first] if pin_first in df.columns else []
    tail = [pin_last] if pin_last in df.columns else []
    return head + middle + tail


def _stacked_bar(df, colors, title, ylabel, out_path):
    fig, ax = plt.subplots(figsize=(13, 7))
    x = np.arange(len(df.index))
    bottom = np.zeros(len(df))
    for col in df.columns:
        vals = df[col].values
        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=colors[col],
            label=col,
            edgecolor="white",
            linewidth=0.3,
        )
        bottom += vals
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12, pad=8)
    ax.margins(x=0.01)
    ax.legend(
        title="Metabolite",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=7,
        title_fontsize=8.5,
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


if __name__ == "__main__":
    main()

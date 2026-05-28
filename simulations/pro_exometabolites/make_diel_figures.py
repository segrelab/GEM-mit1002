#!/usr/bin/env python3
"""
Figures for the diel FBA analysis of Alteromonas MIT1002 growth on Pro exometabolites.

Reads outputs from run_diel_fba.py and produces three figures:

Figure 1 — Headline: predicted MIT1002 growth rate over the diel cycle.
Figure 2 — Substrate hierarchy: which metabolites are binding constraints at each interval.
Figure 3 — Na⁺-cycling activity: flux through the Na⁺-translocating NADH:ubiquinone
           reductase (ec7211_c0) over the diel cycle.

Dark periods: t = 10–22 h and t = 34–46 h (Katie's Figure 4.3 shading).
"""

from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ── Paths ───────────────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).parent
REPO_DIR = SCRIPT_DIR.parent.parent
RESULTS_DIR = SCRIPT_DIR / "results"
CELL_DENSITY_FILE = SCRIPT_DIR / "data/extrapolated_cellcounts.csv"
FIG_DIR = SCRIPT_DIR / "figs"
FIG_DIR.mkdir(exist_ok=True)

# ── Layout constants ────────────────────────────────────────────────────────────

DARK_PERIODS = [(10, 22), (34, 46)]

PALETTE = {
    0.0: "#cccccc",
    0.1: "#a6cee3",
    0.5: "#1f78b4",
    1.0: "#33a02c",
    2.0: "#fb9a99",
    5.0: "#e31a1c",
    10.0: "#ff7f00",
}

# Metabolite display names (short)
MET_LABELS = {
    "EX_cpd00023_e0": "Glu",
    "EX_cpd00080_e0": "G3P",
    "EX_cpd00024_e0": "αKG",
    "EX_cpd00018_e0": "AMP",
    "EX_cpd00091_e0": "UMP",
    "EX_cpd00061_e0": "PEP",
    "EX_cpd00169_e0": "3-PG",
    "EX_cpd00395_e0": "Cys",
    "EX_cpd00477_e0": "NAG",
    "EX_cpd01311_e0": "Dst",
    "EX_cpd00111_e0": "GSSG",
    "EX_cpd02711_e0": "KDPG",
}

NANQR_RXN = "ec7211_c0"


def shade_dark(ax, alpha: float = 0.15) -> None:
    """Shade dark periods on an axis."""
    for d_start, d_end in DARK_PERIODS:
        ax.axvspan(d_start, d_end, color="gray", alpha=alpha, zorder=0, label="_dark")


def load_cell_density() -> pd.DataFrame:
    """Load Pro cell density data (one row per timepoint)."""
    df = pd.read_csv(CELL_DENSITY_FILE)
    return df.drop_duplicates("time_h").sort_values("time_h").reset_index(drop=True)


# ── Figure 1: Headline growth rate ──────────────────────────────────────────────
def figure1(
    growth_df: pd.DataFrame, cell_density: pd.DataFrame, out_path: Path
) -> None:
    """Predicted MIT1002 growth rate vs. time, one line per f, Pro density overlay."""
    fig, ax1 = plt.subplots(figsize=(9, 4))
    ax2 = ax1.twinx()

    # Interval midpoints
    growth_df = growth_df.copy()
    growth_df["midpoint_h"] = (
        growth_df["interval_start_h"] + growth_df["interval_end_h"]
    ) / 2

    # Growth rate lines (one per f)
    f_values = sorted(growth_df["f"].unique())
    for f in f_values:
        sub = growth_df[growth_df["f"] == f].sort_values("midpoint_h")
        mu = sub["growth_rate"].where(sub["status"] == "optimal", other=np.nan)
        ax1.plot(
            sub["midpoint_h"],
            mu,
            "o-",
            color=PALETTE.get(f, "#888888"),
            linewidth=1.5,
            markersize=4,
            label=f"f = {f}",
            zorder=3,
        )

    # Pro cell density on secondary axis
    cd = cell_density[cell_density["time_h"].between(0, 46)]
    ax2.plot(
        cd["time_h"],
        cd["cell_count_mean"] / 1e6,
        "--",
        color="#7570b3",
        linewidth=1.2,
        alpha=0.6,
        label="Pro density",
        zorder=2,
    )
    ax2.fill_between(
        cd["time_h"],
        (cd["cell_count_mean"] - cd["cell_count_sd"].fillna(0)) / 1e6,
        (cd["cell_count_mean"] + cd["cell_count_sd"].fillna(0)) / 1e6,
        color="#7570b3",
        alpha=0.1,
        zorder=1,
    )

    shade_dark(ax1)

    ax1.set_xlabel("Time (h)", fontsize=11)
    ax1.set_ylabel("Predicted growth rate μ (h⁻¹)", fontsize=11)
    ax2.set_ylabel("Pro cell density (×10⁶ cells mL⁻¹)", fontsize=10, color="#7570b3")
    ax2.tick_params(axis="y", labelcolor="#7570b3")
    ax1.set_xlim(0, 46)
    ax1.set_xticks(range(0, 47, 4))

    # Legend: f lines + density + dark period
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles1 = [h for h, l in zip(handles1, labels1) if not l.startswith("_")]
    labels1 = [l for l in labels1 if not l.startswith("_")]
    handles2, labels2 = ax2.get_legend_handles_labels()
    dark_patch = plt.Rectangle((0, 0), 1, 1, fc="gray", alpha=0.3, label="Dark period")
    # Create a figure-level vertical legend placed in the left margin

    ax1.set_title(
        "Predicted Alteromonas MIT1002 growth on Pro-released exometabolites",
        fontsize=11,
        pad=8,
    )
    fig.tight_layout()
    ncol = max(1, len(labels1 + labels2 + ["Dark period"]))
    ncol = min(ncol, 6)
    fig.subplots_adjust(left=0.28)
    fig.legend(
        handles=handles1 + handles2 + [dark_patch],
        labels=labels1 + labels2 + ["Dark period"],
        loc="center left",
        bbox_to_anchor=(0.02, 0.5),
        ncol=1,
        fontsize=8,
        framealpha=0.8,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


# ── Figure 2: Substrate hierarchy ───────────────────────────────────────────────
def figure2(binding_df: pd.DataFrame, out_path: Path, f_rep: float = 1.0) -> None:
    """Heatmap: which exchange bounds are binding at each interval (at f=f_rep)."""
    sub = binding_df[binding_df["f"] == f_rep].copy()
    if sub.empty:
        print(f"No data for f={f_rep} in binding_constraints — skipping Figure 2")
        return

    sub["midpoint_h"] = (sub["interval_start_h"] + sub["interval_end_h"]) / 2
    sub["met_label"] = (
        sub["exchange_reaction"].map(MET_LABELS).fillna(sub["exchange_reaction"])
    )

    # Pivot: rows = metabolite, cols = timepoint, values = |flux|
    piv = sub.pivot_table(
        index="met_label",
        columns="midpoint_h",
        values="flux_at_optimum",
        aggfunc="first",
        fill_value=0,
    )
    piv = piv.abs()  # uptake fluxes are negative; use absolute value

    # Sort rows by mean flux (most-used first)
    piv = piv.loc[piv.mean(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(11, max(3, len(piv) * 0.5 + 1.5)))
    im = ax.imshow(
        piv.values,
        aspect="auto",
        cmap="YlOrRd",
        interpolation="none",
        vmin=0,
    )
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels([f"{t:g}" for t in piv.columns], fontsize=8)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels(piv.index, fontsize=9)
    ax.set_xlabel("Interval midpoint (h)", fontsize=10)
    ax.set_ylabel("Metabolite", fontsize=10)
    ax.set_title(f"Uptake flux magnitude at f = {f_rep} (mmol gDW⁻¹ hr⁻¹)", fontsize=10)

    # Shade dark columns
    col_times = list(piv.columns)
    for d_start, d_end in DARK_PERIODS:
        for ci, t in enumerate(col_times):
            if d_start <= t <= d_end:
                ax.axvline(ci, color="gray", linewidth=4, alpha=0.25, zorder=0)

    plt.colorbar(im, ax=ax, shrink=0.6, label="Uptake flux (mmol gDW⁻¹ hr⁻¹)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


# ── Figure 3: Na⁺ pump activity ─────────────────────────────────────────────────
def figure3(flux_df: pd.DataFrame, out_path: Path) -> None:
    """Flux through Na⁺-translocating NQR (ec7211_c0) vs. time, one line per f."""
    if flux_df.empty:
        print("No flux data — skipping Figure 3")
        return

    nqr_df = flux_df[flux_df["reaction_id"] == NANQR_RXN].copy()
    if nqr_df.empty:
        print(f"Reaction {NANQR_RXN} not found in flux data — skipping Figure 3")
        return

    nqr_df["midpoint_h"] = (nqr_df["interval_start_h"] + nqr_df["interval_end_h"]) / 2

    fig, ax = plt.subplots(figsize=(9, 4))
    shade_dark(ax)

    f_values = sorted(nqr_df["f"].unique())
    for f in f_values:
        sub = nqr_df[nqr_df["f"] == f].sort_values("midpoint_h")
        ax.plot(
            sub["midpoint_h"],
            sub["flux"],
            "o-",
            color=PALETTE.get(f, "#888888"),
            linewidth=1.5,
            markersize=4,
            label=f"f = {f}",
            zorder=3,
        )

    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xlabel("Time (h)", fontsize=11)
    ax.set_ylabel(
        "ec7211_c0 flux (mmol gDW⁻¹ hr⁻¹)\n(NaNQR, Na⁺-translocating NQR)", fontsize=10
    )
    ax.set_xlim(0, 46)
    ax.set_xticks(range(0, 47, 4))

    dark_patch = plt.Rectangle((0, 0), 1, 1, fc="gray", alpha=0.3, label="Dark period")
    handles, labels = ax.get_legend_handles_labels()
    handles = [h for h, l in zip(handles, labels) if not l.startswith("_")]
    labels = [l for l in labels if not l.startswith("_")]
    # Create a figure-level vertical legend placed in the left margin
    ncol = max(1, len(labels + ["Dark period"]))
    ncol = min(ncol, 6)
    fig.tight_layout()
    fig.subplots_adjust(left=0.28)
    fig.legend(
        handles=handles + [dark_patch],
        labels=labels + ["Dark period"],
        loc="center left",
        bbox_to_anchor=(0.02, 0.5),
        ncol=1,
        fontsize=8,
        framealpha=0.8,
    )
    ax.set_title(
        "Na⁺-translocating NQR (ec7211_c0) flux over the diel cycle", fontsize=11, pad=8
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


# Line graph on exchange fluxes over time, one graph per f value, with dark periods shaded
def figure4(flux_df: pd.DataFrame, model: cobra.Model, out_dir: Path) -> None:
    """Line graph: fluxes of all exchange reactions over time, one graph per f."""
    if flux_df.empty:
        print("No flux data — skipping Figure 4")
        return

    # Filter for exchange reactions
    exch_df = flux_df[flux_df["reaction_id"].str.startswith("EX_")].copy()
    if exch_df.empty:
        print("No exchange reactions found in flux data — skipping Figure 4")
        return
    # Add columns for the midpoint hr and the human-readable metabolite label
    exch_df["midpoint_h"] = (
        exch_df["interval_start_h"] + exch_df["interval_end_h"]
    ) / 2
    exch_df["met_name"] = exch_df["reaction_id"].map(
        lambda rxn_id: model.metabolites.get_by_id(rxn_id[3:-3] + "_c0").name
    )

    f_values = sorted(exch_df["f"].unique())
    for f in f_values:
        # Make an outpath for this f value (e.g. "fig4_exchange_fluxes_f1.0.png")
        out_path = out_dir / f"fig4_exchange_fluxes_f{f}.png"
        # Filter for this f value
        sub = exch_df[exch_df["f"] == f].copy()

        # Keep all rows for metabolites that exceed the threshold at any timepoint
        above_thresh = sub.groupby("met_name")["flux"].transform(
            lambda s: s.abs().max() > 0.01
        )
        sub = sub[above_thresh]

        # Make a figure
        fig, ax = plt.subplots(figsize=(9, 4))
        # Shade dark periods
        shade_dark(ax)

        # Plot line for each metabolite
        for met in sub["met_name"].unique():
            met_sub = sub[sub["met_name"] == met].sort_values("midpoint_h")
            ax.plot(
                met_sub["midpoint_h"],
                met_sub["flux"],
                "o-",
                linewidth=1.5,
                markersize=4,
                label=met,
                zorder=3,
            )

        # Style
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xlabel("Time (h)", fontsize=11)
        ax.set_ylabel("Exchange flux (mmol gDW⁻¹ hr⁻¹)", fontsize=10)
        ax.set_xlim(0, 46)
        ax.set_xticks(range(0, 47, 4))

        # Add legend with one entry per metabolite + dark period, placed in the left margin
        dark_patch = plt.Rectangle(
            (0, 0), 1, 1, fc="gray", alpha=0.3, label="Dark period"
        )
        handles, labels = ax.get_legend_handles_labels()
        handles = [h for h, l in zip(handles, labels) if not l.startswith("_")]
        labels = [l for l in labels if not l.startswith("_")]
        ncol = max(1, len(labels + ["Dark period"]))
        ncol = min(ncol, 6)
        fig.tight_layout()
        fig.subplots_adjust(left=0.28)
        fig.legend(
            handles=handles + [dark_patch],
            labels=labels + ["Dark period"],
            loc="center left",
            bbox_to_anchor=(0.02, 0.5),
            ncol=1,
            fontsize=8,
            framealpha=0.8,
        )

        # Save the figure
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved {out_path}")


# ── Main ────────────────────────────────────────────────────────────────────────
def main() -> None:
    growth_file = RESULTS_DIR / "growth_rates.csv"
    flux_file = RESULTS_DIR / "fluxes_long.csv"
    binding_file = RESULTS_DIR / "binding_constraints.csv"

    print("Loading simulation outputs...")
    growth_df = pd.read_csv(growth_file)
    binding_df = pd.read_csv(binding_file)
    flux_df = pd.read_csv(flux_file) if flux_file.exists() else pd.DataFrame()

    # Load the model
    model = cobra.io.read_sbml_model(REPO_DIR / "model.xml")

    cell_density = load_cell_density()

    print(f"  growth_rates: {len(growth_df)} rows")
    print(f"  binding_constraints: {len(binding_df)} rows")
    print(f"  fluxes_long: {len(flux_df)} rows")

    print("\nGenerating Figure 1 (headline growth rates)...")
    figure1(growth_df, cell_density, FIG_DIR / "fig1_growth_rates.png")

    print("Generating Figure 2 (substrate hierarchy at f=1.0)...")
    figure2(binding_df, FIG_DIR / "fig2_substrate_hierarchy.png", f_rep=1.0)

    print("Generating Figure 3 (NaNQR flux)...")
    figure3(flux_df, FIG_DIR / "fig3_nanqr_flux.png")

    print("Generating Figure 4 (exchange fluxes over time)...")
    figure4(flux_df, model, FIG_DIR)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Integrated diel dynamics for Pro/Amac co-culture (f = 10, NH3 removed from BASAL).

Plots five lines on a single panel with three y-axes:

  1. Glutamate concentration — experimental (smoothed, with replicate scatter)
  2. Glutamate concentration — simulated (experimental − cumulative Amac uptake)
  3. Ammonium concentration — predicted (integrated Amac NH3 release; starts at 0)
  4. Prochlorococcus cell density (input data, ×10⁶ cells mL⁻¹)
  5. Alteromonas biomass (μg DW L⁻¹), integrated from FBA growth rates

Integration details
-------------------
Amac biomass evolves as dX/dt = μ(t) · X(t). With μ piecewise-constant on each
2-h interval:
    X(t_{i+1}) = X(t_i) · exp(μ_i · Δt)

Time-integrated biomass over an interval (needed to convert per-gDW fluxes to
per-L mass changes):
    ∫ X(s) ds  =  X(t_i) · [exp(μ_i Δt) − 1] / μ_i   (μ_i > 0)
                 X(t_i) · Δt                          (μ_i ≈ 0)

Per-interval glutamate uptake (mmol/L) = (−v_glu) · ∫X ds, where v_glu is the
EX_cpd00023_e0 flux from FBA (negative = uptake under COBRA convention).
Cumulative uptake is converted to nM (×1e6) and subtracted from the
experimental concentration to get the predicted with-Amac concentration. As
Helen noted, this can go negative when Pro starts reabsorbing glutamate at
night — that's expected and worth showing.

NH3 release uses the same integration but with v_nh3 directly (positive flux =
secretion). Starts from 0 nM (BASAL has NH3 removed in the simulation).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Paths and parameters ────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results"
DATA_DIR = SCRIPT_DIR / "data"
FIG_DIR = SCRIPT_DIR / "figs"
FIG_DIR.mkdir(exist_ok=True)

F_PLOT = 10.0
ALT_DW_G = 2.5e-13  # g/cell
# Initial Amac biomass is set to be self-consistent with f: the FBA bound
# at f=10 implicitly assumes 10 Pro per Amac, so Amac inoculum = Pro(t=0) / f.
# Computed below from the measured initial Pro density.

DT_H = 2.0
DARK_PERIODS = [(10, 22), (34, 46)]
SMOOTHING_WINDOW = 3

GLU_EX = "EX_cpd00023_e0"
NH3_EX = "EX_cpd00013_e0"


def shade_dark(ax, alpha: float = 0.15) -> None:
    for d0, d1 in DARK_PERIODS:
        ax.axvspan(d0, d1, color="gray", alpha=alpha, zorder=0)


def main() -> None:
    growth_df = pd.read_csv(RESULTS_DIR / "growth_rates.csv")
    fluxes_df = pd.read_csv(RESULTS_DIR / "fluxes_long.csv")
    glu_means = pd.read_csv(DATA_DIR / "ProDiel_filtered_meanByTimepoint.csv")
    pro_density = pd.read_csv(DATA_DIR / "extrapolated_cellcounts.csv")
    glu_reps_path = DATA_DIR / "ProDiel_filtered_replicates.csv"
    glu_reps = pd.read_csv(glu_reps_path) if glu_reps_path.exists() else None

    # f=10 subsets, sorted by interval start
    g = (
        growth_df[growth_df["f"] == F_PLOT]
        .sort_values("interval_start_h")
        .reset_index(drop=True)
    )
    fx = fluxes_df[fluxes_df["f"] == F_PLOT]
    v_glu_by_t = fx[fx["reaction_id"] == GLU_EX].set_index("interval_start_h")["flux"]
    v_nh3_by_t = fx[fx["reaction_id"] == NH3_EX].set_index("interval_start_h")["flux"]

    # Time grid: interval starts plus the final end
    t_starts = g["interval_start_h"].to_numpy()
    t_ends = g["interval_end_h"].to_numpy()
    n = len(t_starts)
    times = np.concatenate(([t_starts[0]], t_ends))  # length n+1

    mu = g["growth_rate"].fillna(0.0).to_numpy()

    # Self-consistent initial Amac biomass: Pro(t=0) / f cells/mL
    pro_dedup = pro_density.drop_duplicates("time_h").sort_values("time_h")
    init_pro_cells_per_ml = float(pro_dedup["cell_count_mean"].iloc[0])
    init_amac_cells_per_ml = init_pro_cells_per_ml / F_PLOT
    X0_GDW_PER_L = init_amac_cells_per_ml * ALT_DW_G * 1000

    # Amac biomass at each timepoint
    X = np.empty(n + 1)
    X[0] = X0_GDW_PER_L
    for i in range(n):
        X[i + 1] = X[i] * np.exp(mu[i] * DT_H)

    # Time-integrated biomass over each interval
    int_X = np.empty(n)
    for i in range(n):
        if mu[i] > 1e-12:
            int_X[i] = X[i] * (np.exp(mu[i] * DT_H) - 1.0) / mu[i]
        else:
            int_X[i] = X[i] * DT_H

    # Per-interval fluxes; missing (e.g. infeasible) treated as zero
    glu_flux = np.array([v_glu_by_t.get(t, 0.0) for t in t_starts])
    nh3_flux = np.array([v_nh3_by_t.get(t, 0.0) for t in t_starts])

    # Cumulative glutamate uptake (mmol/L → nM)
    glu_uptake_per_int = -glu_flux * int_X  # mmol/L (uptake is positive here)
    glu_uptake_cum_nM = np.concatenate(([0.0], np.cumsum(glu_uptake_per_int))) * 1e6

    # Experimental glutamate (smooth here so we don't depend on convert_data_to_rates outputs)
    glu_exp = (
        glu_means[glu_means["CleanName"] == "glutamic_acid"]
        .sort_values("timepoint")
        .reset_index(drop=True)
    )
    glu_exp["smoothed_nM"] = (
        glu_exp["mean_nM"].rolling(SMOOTHING_WINDOW, center=True, min_periods=1).mean()
    )
    smooth_lookup = glu_exp.set_index("timepoint")["smoothed_nM"]
    glu_exp_at_times = np.array([smooth_lookup.get(t, np.nan) for t in times])

    # Simulated glutamate = experimental − cumulative Amac uptake
    glu_sim = glu_exp_at_times - glu_uptake_cum_nM

    # Amac biomass in μg/L
    X_ugL = X * 1e6

    # Pro density (one row per timepoint)
    pro = (
        pro_density.drop_duplicates("time_h")
        .sort_values("time_h")
        .reset_index(drop=True)
    )

    # Diagnostic print
    print(f"Initial Pro density:         {init_pro_cells_per_ml:.2e} cells/mL")
    print(f"Initial Amac density (= Pro/f): {init_amac_cells_per_ml:.2e} cells/mL")
    print(
        f"Initial Amac biomass:        {X0_GDW_PER_L*1e6:.2f} μg/L "
        f"({ALT_DW_G*1e15:.0f} fg/cell)"
    )
    print(f"Final Amac biomass:          {X[-1]*1e6:.2f} μg/L")
    print(f"Total glutamate consumed:    {glu_uptake_cum_nM[-1]:.1f} nM")
    print(
        f"Peak experimental glutamate: {np.nanmax(glu_exp_at_times):.1f} nM "
        f"at t={times[int(np.nanargmax(glu_exp_at_times))]:g} h"
    )

    # ── Plot ────────────────────────────────────────────────────────────────────
    fig, ax1 = plt.subplots(figsize=(11.5, 5.2))
    ax2 = ax1.twinx()  # Pro density (×10⁶ cells/mL)
    ax3 = ax1.twinx()  # Amac biomass (μg/L)
    ax3.spines["right"].set_position(("outward", 60))

    shade_dark(ax1)
    ax1.axhline(0, color="black", lw=0.5, zorder=1)

    c_glu_exp = "#1f78b4"  # blue (experimental)
    c_glu_sim = "#a6cee3"  # light blue (simulated)
    c_nh4 = "#e31a1c"  # red (NH4)
    c_pro = "#7570b3"  # purple (Pro)
    c_amac = "#33a02c"  # green (Amac)

    ax1.plot(
        times,
        glu_exp_at_times,
        "o-",
        color=c_glu_exp,
        lw=2,
        ms=5,
        label="Glu — experimental (smoothed)",
        zorder=4,
    )
    ax1.plot(
        times,
        glu_sim,
        "s--",
        color=c_glu_sim,
        lw=2,
        ms=4,
        label="Glu — simulated (exp − Amac uptake)",
        zorder=4,
    )

    ax1.set_xlabel("Time (h)", fontsize=11)
    ax1.set_ylabel("Concentration (nM)", fontsize=11)
    ax1.set_xlim(0, 46)
    ax1.set_xticks(range(0, 47, 4))

    # Pro density (right inner)
    ax2.plot(
        pro["time_h"],
        pro["cell_count_mean"] / 1e6,
        ":",
        color=c_pro,
        lw=2.2,
        label="Pro density",
    )
    if "cell_count_sd" in pro.columns:
        ax2.fill_between(
            pro["time_h"],
            (pro["cell_count_mean"] - pro["cell_count_sd"].fillna(0)) / 1e6,
            (pro["cell_count_mean"] + pro["cell_count_sd"].fillna(0)) / 1e6,
            color=c_pro,
            alpha=0.10,
        )
    ax2.set_ylabel("Pro density (×10⁶ cells mL⁻¹)", color=c_pro, fontsize=10)
    ax2.tick_params(axis="y", labelcolor=c_pro)

    # Amac biomass (right outer)
    ax3.plot(times, X_ugL, "-", color=c_amac, lw=2.6, label="Amac biomass")
    ax3.set_ylabel("Amac biomass (μg DW L⁻¹)", color=c_amac, fontsize=10)
    ax3.tick_params(axis="y", labelcolor=c_amac)

    # Combined legend
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    h3, l3 = ax3.get_legend_handles_labels()
    dark_patch = plt.Rectangle((0, 0), 1, 1, fc="gray", alpha=0.3, label="Dark period")
    ax1.legend(
        handles=h1 + h2 + h3 + [dark_patch],
        labels=l1 + l2 + l3 + ["Dark period"],
        fontsize=8,
        loc="upper right",
        ncol=2,
        framealpha=0.85,
    )

    fig.suptitle(
        f"Diel dynamics: glutamate, ammonium, Pro, Amac (f = {F_PLOT:g}, NH₃ removed from medium)",
        fontsize=11,
    )
    plt.tight_layout()
    out = FIG_DIR / "fig_integrated_dynamics.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()

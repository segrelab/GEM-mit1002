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
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Paths and parameters ────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS_DIR = SCRIPT_DIR / "results"
DATA_DIR = SCRIPT_DIR / "data"
FIG_DIR = SCRIPT_DIR / "figs"
FIG_DIR.mkdir(exist_ok=True)

# Import the shared plot styles from tools/
sys.path.append(str(REPO_ROOT))
from tools.plot_styles import set_plot_style, summer_colors

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
        ax.axvspan(d0, d1, color=summer_colors["dark_tan"], alpha=alpha, zorder=0)


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

    # Cumulative ammonium release (mmol/L → nM), starting from 0
    nh3_release_per_int = nh3_flux * int_X
    nh4_cum_nM = np.concatenate(([0.0], np.cumsum(nh3_release_per_int))) * 1e6

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
    # Two stacked panels sharing the time axis:
    #   top    — experimental data only (glutamate + Pro cell density)
    #   bottom — simulation outputs (glutamate + ammonium + Amac biomass)
    fig, (ax_exp, ax_sim) = plt.subplots(
        2, 1, figsize=(7.2, 7.2), sharex=True, gridspec_kw={"hspace": 0.18}
    )
    ax_exp_r = ax_exp.twinx()  # Pro density (×10⁶ cells/mL)
    ax_sim_r = ax_sim.twinx()  # Amac biomass (μg/L)

    # Background shading for dark periods
    for ax in (ax_exp, ax_sim):
        shade_dark(ax, alpha=0.4)
        ax.axhline(0, color="gray", lw=0.5, zorder=1)

    # Defnine colors for each line
    c_glu_exp = summer_colors["dark_pink"]
    c_glu_sim = summer_colors["pink"]
    c_nh4 = summer_colors["yellow"]
    c_pro = summer_colors["green"]
    c_amac = summer_colors["teal"]

    # Only the first diel cycle is shown; y-limits are set from this window so
    # the (unplotted) second day can't inflate them
    x_max = 22
    in_win = times <= x_max

    # ── Top panel: experimental data ────────────────────────────────────────────
    ax_exp.plot(
        times,
        glu_exp_at_times,
        "o-",
        color=c_glu_exp,
        lw=2,
        ms=5,
        label="Glutamate (Experimental)",
        zorder=4,
    )
    ax_exp.set_ylabel("Concentration (nM)", fontsize=11)
    # Headroom above the data so the upper-left legend never sits on a line
    ax_exp.set_ylim(0, np.nanmax(glu_exp_at_times[in_win]) * 1.45)

    ax_exp_r.plot(
        pro["time_h"],
        pro["cell_count_mean"] / 1e6,
        ":",
        color=c_pro,
        lw=2.2,
        label="Pro density",
    )
    if "cell_count_sd" in pro.columns:
        ax_exp_r.fill_between(
            pro["time_h"],
            (pro["cell_count_mean"] - pro["cell_count_sd"].fillna(0)) / 1e6,
            (pro["cell_count_mean"] + pro["cell_count_sd"].fillna(0)) / 1e6,
            color=c_pro,
            alpha=0.10,
        )
    # Set the y-axis limits so that the first day's value take up most of the space (the second day is just a repeat of the first)
    ax_exp_r.set_ylim(0, 190)
    ax_exp_r.set_ylabel(
        "Prochlorococcus density\n(×10⁶ cells mL⁻¹)", color=c_pro, fontsize=10
    )
    ax_exp_r.tick_params(axis="y", labelcolor=c_pro)

    # ── Bottom panel: simulation outputs ────────────────────────────────────────
    # Glutamate and ammonium are both nM, so they share the left axis
    ax_sim.plot(
        times,
        glu_sim,
        "s--",
        color=c_glu_sim,
        lw=2,
        ms=4,
        label="Glutamate (Simulated)",
        zorder=4,
    )
    ax_sim.plot(
        times,
        nh4_cum_nM,
        "^-",
        color=c_nh4,
        lw=2,
        ms=4,
        label="Ammonium (Simulated)",
        zorder=4,
    )
    ax_sim.set_ylabel("Concentration (nM)", fontsize=11)
    # Subset the y-axis to only show positive concentrations
    # The simulated glutamate can go negative when Pro starts reabsorbing it
    sim_top = max(np.nanmax(glu_sim[in_win]), np.nanmax(nh4_cum_nM[in_win]))
    ax_sim.set_ylim(0, sim_top * 1.55)

    ax_sim_r.plot(times, X_ugL, "-", color=c_amac, lw=2.6, label="Amac biomass")
    x_lo, x_hi = X_ugL[in_win].min(), X_ugL[in_win].max()
    x_span = x_hi - x_lo
    ax_sim_r.set_ylim(x_lo - 0.08 * x_span, x_hi + 0.55 * x_span)
    ax_sim_r.set_ylabel(
        "A. macleodii biomass\n(μg DW L⁻¹)", color=c_amac, fontsize=10
    )
    ax_sim_r.tick_params(axis="y", labelcolor=c_amac)

    # Shared x-axis: only label the bottom panel
    # Subset the x-axis (time) to the first diel cycle
    ax_sim.set_xlabel("Time (h)", fontsize=11)
    ax_sim.set_xlim(0, x_max)
    ax_sim.set_xticks(range(0, x_max + 1, 4))

    # Per-panel legends, in the headroom cleared above the data
    dark_patch = plt.Rectangle(
        (0, 0), 1, 1, fc=summer_colors["dark_tan"], label="Dark period"
    )
    for ax_l, ax_r in ((ax_exp, ax_exp_r), (ax_sim, ax_sim_r)):
        hl, ll = ax_l.get_legend_handles_labels()
        hr, lr = ax_r.get_legend_handles_labels()
        ax_l.legend(
            handles=hl + hr + [dark_patch],
            labels=ll + lr + ["Dark period"],
            fontsize=8,
            loc="upper left",
            ncol=1,
            framealpha=0.85,
        )

    # Gray axes/ticks/labels. set_plot_style() drops the right spine, which the
    # twinned axes need, so put theirs back (and drop their duplicate left one).
    for ax in (ax_exp, ax_exp_r, ax_sim, ax_sim_r):
        set_plot_style(ax)
    for ax in (ax_exp_r, ax_sim_r):
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color("gray")
        ax.spines["left"].set_visible(False)
    # set_plot_style() grays every label; restore the per-series colors that tie
    # each right-hand axis to its line
    ax_exp_r.yaxis.label.set_color(c_pro)
    ax_sim_r.yaxis.label.set_color(c_amac)

    fig.suptitle(
        f"Diel dynamics: glutamate, ammonium, Pro, Amac (f = {F_PLOT:g}, NH₃ removed from medium)",
        fontsize=11,
        color="gray",
    )
    # tight_layout is skipped: it warns on twinx axes; hspace + bbox_inches handle spacing
    out = FIG_DIR / "fig_integrated_dynamics.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Integrated diel dynamics for Pro/Amac co-culture (f = 10, NH3 removed from BASAL).

Plots two panels:
The top panel shows the experimental data only:
  1. Glutamate concentration — experimental (smoothed)
  2. Prochlorococcus cell density (input data, ×10⁶ cells mL⁻¹)
The bottom panel shows the simulation outputs:
  1. Glutamate concentration — simulated (experimental − cumulative Amac uptake)
  2. Ammonium concentration — predicted (integrated Amac NH3 release; starts at 0)
  3. Alteromonas biomass (μg DW L⁻¹), integrated from FBA growth rates

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

import sys
from pathlib import Path

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

# Figure geometry, in inches. FIG_W_IN and PANEL_H_IN are the *total* footprint
# of each panel — axis labels and the right-hand legend included — so the figure
# is FIG_W_IN wide by 2 * PANEL_H_IN tall. Everything else (the margins, the
# width of the legend column) is measured from the rendered figure by
# fit_layout() below, so the panels stay as wide as the labels allow and no
# whitespace is left over when the labels change length.
FIG_W_IN = 5.36
PANEL_H_IN = 1.7
PANEL_GAP_IN = 0.20
LEGEND_GAP_IN = 0.10  # between the right-hand y-axis label and the legend
PAD_IN = 0.04  # outer padding on the other three sides

GLU_EX = "EX_cpd00023_e0"
NH3_EX = "EX_cpd00013_e0"


def shade_dark(ax, alpha: float = 0.15) -> None:
    for d0, d1 in DARK_PERIODS:
        ax.axvspan(d0, d1, color=summer_colors["dark_tan"], alpha=alpha, zorder=0)


def fit_layout(fig, panels) -> None:
    """Size the axes so the figure is exactly filled, whatever the labels say.

    `panels` is a list of (left_axis, right_axis, legend). Each item's decoration
    widths (tick labels, y-axis labels, legend) are measured from a rendered
    figure, and the axes are then stretched to take up whatever is left. Called
    a few times because resizing the axes changes the tick labels, which changes
    the measurements; it settles after two or three passes.
    """
    fig_w, fig_h = fig.get_size_inches()
    for _ in range(3):
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()

        def span(artist, tight=True):
            """Bounding box of an artist, in inches: (x0, x1, y0, y1)."""
            # Axis objects need get_tightbbox() — get_window_extent() leaves out
            # the tick labels and the axis label, which is exactly what we're
            # trying to measure here.
            bb = artist.get_tightbbox(rend) if tight else None
            if bb is None:  # nothing drawn (e.g. a shared, label-less x-axis)
                bb = artist.get_window_extent(rend)
            return bb.x0 / fig.dpi, bb.x1 / fig.dpi, bb.y0 / fig.dpi, bb.y1 / fig.dpi

        left_axes = [ax_l for ax_l, _, _ in panels]
        pos = [ax.get_position() for ax in left_axes]
        ax_x0 = min(pp.x0 for pp in pos) * fig_w
        ax_x1 = max(pp.x1 for pp in pos) * fig_w
        ax_y0 = min(pp.y0 for pp in pos) * fig_h

        # Width taken by the y-axis labels on either side, and by the legend
        left_deco = max(ax_x0 - span(ax.yaxis)[0] for ax in left_axes)
        right_deco = max(span(ax_r.yaxis)[1] - ax_x1 for _, ax_r, _ in panels)
        legend_w = max(
            span(leg, tight=False)[1] - span(leg, tight=False)[0]
            for _, _, leg in panels
        )
        bottom_deco = ax_y0 - span(left_axes[-1].xaxis)[2]

        axes_w = (
            fig_w - left_deco - right_deco - LEGEND_GAP_IN - legend_w - 2 * PAD_IN
        )
        axes_h = (fig_h - bottom_deco - 2 * PAD_IN - PANEL_GAP_IN) / len(panels)
        fig.subplots_adjust(
            left=(PAD_IN + left_deco) / fig_w,
            right=(PAD_IN + left_deco + axes_w) / fig_w,
            bottom=(PAD_IN + bottom_deco) / fig_h,
            top=(fig_h - PAD_IN) / fig_h,
            hspace=PANEL_GAP_IN / axes_h,
        )

        # Park the legend column just right of the widest y-axis label
        legend_x = (PAD_IN + left_deco + axes_w + right_deco + LEGEND_GAP_IN) / fig_w
        for ax_l, _, leg in panels:
            leg.set_bbox_to_anchor((legend_x, ax_l.get_position().y1), fig.transFigure)


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
    #
    # The panels are wide and short, so the legends sit outside on the right;
    # the margins reserve room for them inside the FIG_W_IN total width.
    fig_w = FIG_W_IN
    fig_h = 2 * PANEL_H_IN
    fig, (ax_exp, ax_sim) = plt.subplots(2, 1, figsize=(fig_w, fig_h), sharex=True)
    ax_exp_r = ax_exp.twinx()  # Pro density (×10⁶ cells/mL)
    ax_sim_r = ax_sim.twinx()  # Amac biomass (μg/L)

    # Background shading for dark periods
    for ax in (ax_exp, ax_sim):
        shade_dark(ax, alpha=0.4)
        ax.axhline(0, color="gray", lw=0.5, zorder=1)

    # Defnine colors for each line
    c_glu = summer_colors["dark_pink"]
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
        color=c_glu,
        lw=1.4,
        ms=3,
        label="Glutamate",
        zorder=4,
    )
    ax_exp.set_ylabel("Concentration (nM)", fontsize=7)
    ax_exp.set_ylim(0, np.nanmax(glu_exp_at_times[in_win]) * 1.12)

    ax_exp_r.plot(
        pro["time_h"],
        pro["cell_count_mean"] / 1e6,
        "-",
        color=c_pro,
        lw=1.6,
        label="Prochlorococcus",
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
    ax_exp_r.set_ylim(0, 130)
    ax_exp_r.set_ylabel("Prochlorococcus\n(×10⁶ cells mL⁻¹)", color=c_pro, fontsize=7)
    ax_exp_r.tick_params(axis="y", labelcolor=c_pro)

    # ── Bottom panel: simulation outputs ────────────────────────────────────────
    # Glutamate and ammonium are both nM, so they share the left axis
    ax_sim.plot(
        times,
        glu_sim,
        "o--",
        color=c_glu,
        lw=1.4,
        ms=2.8,
        label="Glutamate",
        zorder=4,
    )
    ax_sim.plot(
        times,
        nh4_cum_nM,
        "s--",
        color=c_nh4,
        lw=1.4,
        ms=2.8,
        label="Ammonium",
        zorder=4,
    )
    ax_sim.set_ylabel("Concentration (nM)", fontsize=7)
    # Subset the y-axis to only show positive concentrations
    # The simulated glutamate can go negative when Pro starts reabsorbing it
    sim_top = max(np.nanmax(glu_sim[in_win]), np.nanmax(nh4_cum_nM[in_win]))
    ax_sim.set_ylim(0, sim_top * 1.12)

    ax_sim_r.plot(times, X_ugL, "--", color=c_amac, lw=1.8, label="MIT1002")
    x_lo, x_hi = X_ugL[in_win].min(), X_ugL[in_win].max()
    x_span = x_hi - x_lo
    ax_sim_r.set_ylim(x_lo - 0.08 * x_span, x_hi + 0.12 * x_span)
    ax_sim_r.set_ylabel("MIT1002\n(μg DW L⁻¹)", color=c_amac, fontsize=7)
    ax_sim_r.tick_params(axis="y", labelcolor=c_amac)

    # Shared x-axis: only label the bottom panel
    # Subset the x-axis (time) to the first diel cycle
    ax_sim.set_xlabel("Time (h)", fontsize=7)
    ax_sim.set_xlim(0, x_max)
    ax_sim.set_xticks(range(0, x_max + 1, 4))
    for ax in (ax_exp, ax_exp_r, ax_sim, ax_sim_r):
        ax.tick_params(labelsize=6, length=2, pad=1.5)

    # Per-panel legends. fit_layout() places them in a column to the right of
    # the panels, anchored in figure coordinates so the two columns line up.
    dark_patch = plt.Rectangle(
        (0, 0), 1, 1, fc=summer_colors["dark_tan"], label="Dark period"
    )
    panels = []
    for ax_l, ax_r in ((ax_exp, ax_exp_r), (ax_sim, ax_sim_r)):
        hl, ll = ax_l.get_legend_handles_labels()
        hr, lr = ax_r.get_legend_handles_labels()
        leg = ax_l.legend(
            handles=hl + hr + [dark_patch],
            labels=ll + lr + ["Dark period"],
            fontsize=6,
            loc="upper left",
            bbox_transform=fig.transFigure,
            borderaxespad=0,
            handlelength=1.6,
            handletextpad=0.5,
            labelspacing=0.4,
            ncol=1,
            frameon=False,
        )
        panels.append((ax_l, ax_r, leg))

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

    # Measure the labels and legends, then stretch the axes to fill the figure
    fit_layout(fig, panels)

    out = FIG_DIR / "fig_integrated_dynamics.png"
    # Note: no bbox_inches="tight" — it would re-crop away the fixed figure size
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Figure: {fig_w:.2f} x {fig_h:.2f} in ({PANEL_H_IN:.2f} in per panel)")
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()

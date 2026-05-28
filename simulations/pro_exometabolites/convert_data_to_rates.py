"""
Compute per-Prochlorococcus-cell metabolite release rates from a diel timecourse.

For each metabolite, smooths the concentration timecourse with a centered
rolling mean (default 3-point window) to suppress measurement noise, then
computes dC/dt between consecutive timepoints (forward differences). The
dC/dt is normalized by the mean Pro cell density across the interval and
negative values are clipped to zero. Clipping treats intervals where Pro is
net-reabsorbing a metabolite as supplying zero flux to a heterotrophic
neighbor (a conservative lower bound).

Smoothing is necessary because the concentration measurements have ~30% RSD
(see Katie's Table 4.1), and differencing noisy data amplifies the noise:
the noise in C(t+1) - C(t) is sqrt(2) times the noise in either point alone.
Without smoothing, the dC/dt timecourse oscillates between release and
reabsorption on a 2-hour timescale, which is faster than the underlying
biology and reflects measurement noise rather than real Pro physiology.

The output rate at row i is the rate for the interval [t_i, t_{i+1}], assigned
to the start of the interval. Units are nmol per Pro cell per hour. This is
NOT yet normalized to Alteromonas dry weight; that conversion happens later
when these rates are turned into FBA exchange bounds.

Inputs:
    ProDiel_filtered_meanByTimepoint.csv: filtered exometabolite concentrations
        (long format: CleanName, timepoint, mean_nM, sd_nM, n_reps)
    extrapolated_cellcounts.csv: Pro cell density timecourse from Natalie
    ProDiel_filtered_replicates.csv (optional, for diagnostic plot only):
        individual replicate measurements, used to show scatter on plot

Output:
    ProDiel_per_pro_cell_rates.csv: long format with columns
        metabolite, interval_start_h, interval_end_h, dC_dt_nM_per_hr,
        mean_pro_density_cells_per_mL, per_pro_cell_rate_nmol_per_hr,
        per_pro_cell_rate_nmol_per_hr_raw, cell_density_quality
    ProDiel_per_pro_cell_rates.png: diagnostic plot showing for each metabolite
        the raw and smoothed concentrations (top) and rates computed from both
        (bottom). Lets you visually confirm smoothing is reasonable.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# Constants for file paths
FILE_DIR = Path(__file__).parent
OUT_DIR = Path(FILE_DIR, "results")
# If the output directory doesn't exist, create it
OUT_DIR.mkdir(exist_ok=True)

# Set the paths/names for the input and output files
INPUT_METABOLITES = Path(FILE_DIR, "data", "ProDiel_filtered_meanByTimepoint.csv")
INPUT_CELLCOUNTS = Path(FILE_DIR, "data", "extrapolated_cellcounts.csv")
INPUT_REPLICATES = Path(FILE_DIR, "data", "ProDiel_filtered_replicates.csv")
OUTPUT_FILE = Path(OUT_DIR, "ProDiel_per_pro_cell_rates.csv")
PLOT_FILE = Path(OUT_DIR, "ProDiel_per_pro_cell_rates.png")

# Define the experiments dark periods (from Katie's chapter Figure 4.3 shading)
DARK_PERIODS = [(10, 22), (34, 46)]

# Smoothing window (number of timepoints) for centered rolling mean.
# 3 = mild smoothing, preserves diel features. Set to 1 to disable smoothing.
SMOOTHING_WINDOW = 3


def main() -> None:
    # Load the data
    metabolites = load_metabolites(INPUT_METABOLITES)
    cell_densities = load_cell_densities(INPUT_CELLCOUNTS)

    # Add a smoothed-concentration column for use in dC/dt and the plot
    metabolites = smooth_concentrations(metabolites, window=SMOOTHING_WINDOW)

    # Calculate per-Pro-cell rates from smoothed concentrations
    rates = compute_per_pro_cell_rates(metabolites, cell_densities)

    # Save the results
    rates.to_csv(OUTPUT_FILE, index=False)
    print(f"Wrote {OUTPUT_FILE}")

    # Plot diagnostic (raw + smoothed concs and rates)
    replicates = load_replicates(INPUT_REPLICATES) if INPUT_REPLICATES.exists() else None
    plot_rates_diagnostic(metabolites, cell_densities, replicates, PLOT_FILE)


def load_cell_densities(path: Path) -> pd.DataFrame:
    """Load cell density data and deduplicate to one row per timepoint.

    cell_count_mean is the same across the three replicates at each timepoint,
    so we keep just one representative row per time.
    """
    df = pd.read_csv(path)
    dedup = (
        df.drop_duplicates("time_h")[
            ["time_h", "cell_count_mean", "cell_count_sd", "typeOfData"]
        ]
        .sort_values("time_h")
        .reset_index(drop=True)
    )
    return dedup


def load_metabolites(path: Path) -> pd.DataFrame:
    """Load filtered metabolite concentrations."""
    df = pd.read_csv(path)
    return df.sort_values(["CleanName", "timepoint"]).reset_index(drop=True)


def load_replicates(path: Path) -> pd.DataFrame:
    """Load individual replicate concentrations (for diagnostic plot only)."""
    return pd.read_csv(path)


def smooth_concentrations(
    metabolites: pd.DataFrame, window: int = SMOOTHING_WINDOW
) -> pd.DataFrame:
    """Add a `mean_nM_smoothed` column via centered rolling mean per metabolite.

    Uses min_periods=1 so endpoints get a smaller-window average rather than NaN.
    With window=1, the smoothed values equal the raw means (smoothing disabled).
    """
    df = metabolites.copy()
    df["mean_nM_smoothed"] = (
        df.sort_values("timepoint")
        .groupby("CleanName")["mean_nM"]
        .transform(
            lambda s: s.rolling(window=window, center=True, min_periods=1).mean()
        )
    )
    return df


def compute_per_pro_cell_rates(
    metabolites: pd.DataFrame,
    cell_densities: pd.DataFrame,
) -> pd.DataFrame:
    """Compute per-Pro-cell release rates via forward differences of smoothed concs.

    For each consecutive pair of timepoints (t_i, t_{i+1}) at which a metabolite
    is measured, computes:
        dC/dt        = (Cs(t_{i+1}) - Cs(t_i)) / (t_{i+1} - t_i)  [nM/hr]
                       where Cs is the smoothed concentration
        N_mean       = (N(t_i) + N(t_{i+1})) / 2                   [cells/mL]
        rate_raw     = dC_dt / (N_mean * 1000)                     [nmol/cell/hr]
                       (the 1000 converts cells/mL -> cells/L)
        rate_clipped = max(rate_raw, 0)
    """
    rows = []
    cd_lookup = cell_densities.set_index("time_h")

    # Column to use for dC/dt; falls back to raw if smoothed missing
    conc_col = "mean_nM_smoothed" if "mean_nM_smoothed" in metabolites.columns else "mean_nM"

    for metabolite, group in metabolites.groupby("CleanName"):
        g = group.sort_values("timepoint").reset_index(drop=True)

        for i in range(len(g) - 1):
            t_start = g.loc[i, "timepoint"]
            t_end = g.loc[i + 1, "timepoint"]
            dt = t_end - t_start

            c_start_nM = g.loc[i, conc_col]
            c_end_nM = g.loc[i + 1, conc_col]
            dC_dt_nM_per_hr = (c_end_nM - c_start_nM) / dt

            # Look up cell densities at both endpoints
            if t_start not in cd_lookup.index or t_end not in cd_lookup.index:
                continue

            n_start = cd_lookup.loc[t_start, "cell_count_mean"]
            n_end = cd_lookup.loc[t_end, "cell_count_mean"]
            type_start = cd_lookup.loc[t_start, "typeOfData"]
            type_end = cd_lookup.loc[t_end, "typeOfData"]

            mean_density_cells_per_mL = (n_start + n_end) / 2
            mean_density_cells_per_L = mean_density_cells_per_mL * 1000

            # nM/hr = nmol/L/hr; divide by cells/L -> nmol/(cell*hr)
            rate_raw = dC_dt_nM_per_hr / mean_density_cells_per_L
            rate_clipped = max(rate_raw, 0)

            quality = (
                "measured"
                if type_start == "measured" and type_end == "measured"
                else "partial_estimate"
            )

            rows.append(
                {
                    "metabolite": metabolite,
                    "interval_start_h": t_start,
                    "interval_end_h": t_end,
                    "dC_dt_nM_per_hr": dC_dt_nM_per_hr,
                    "mean_pro_density_cells_per_mL": mean_density_cells_per_mL,
                    "per_pro_cell_rate_nmol_per_hr": rate_clipped,
                    "per_pro_cell_rate_nmol_per_hr_raw": rate_raw,
                    "cell_density_quality": quality,
                }
            )

    return pd.DataFrame(rows)


def _rate_from_concs(
    times: np.ndarray, concs: np.ndarray, density: np.ndarray
) -> np.ndarray:
    """Forward-difference per-Pro-cell rate at each interval.

    Returns array of length len(times) - 1, in nmol/(cell*hr).
    """
    dC_dt = np.diff(concs) / np.diff(times)  # nM/hr = nmol/L/hr
    mean_density_cells_per_L = ((density[:-1] + density[1:]) / 2) * 1000
    return dC_dt / mean_density_cells_per_L


def plot_rates_diagnostic(
    metabolites: pd.DataFrame,
    cell_densities: pd.DataFrame,
    replicates: pd.DataFrame | None,
    output_path: Path,
) -> None:
    """For each metabolite, plot raw + smoothed concentrations and the rates
    computed from each. Lets you see what smoothing is doing and confirm the
    underlying biology is preserved.

    Layout: 4 columns, with each metabolite occupying a 2-row block
    (concentration on top, rate below).
    """
    mets_sorted = sorted(metabolites["CleanName"].unique())
    n_mets = len(mets_sorted)
    ncols = 4
    n_met_rows = (n_mets + ncols - 1) // ncols
    nrows = n_met_rows * 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 2.2 * nrows))

    cd_lookup = cell_densities.set_index("time_h")["cell_count_mean"].to_dict()

    for i, met in enumerate(mets_sorted):
        met_row = i // ncols
        col = i % ncols
        ax_conc = axes[met_row * 2, col]
        ax_rate = axes[met_row * 2 + 1, col]

        # Pull this metabolite's data
        mdata = metabolites[metabolites["CleanName"] == met].sort_values("timepoint")
        times = mdata["timepoint"].values
        raw_means = mdata["mean_nM"].values
        smoothed = mdata["mean_nM_smoothed"].values

        # --- Concentration panel ---
        if replicates is not None:
            rep_data = replicates[replicates["CleanName"] == met]
            ax_conc.scatter(
                rep_data["timepoint"], rep_data["nM"],
                color="#666", alpha=0.5, s=12, zorder=2,
            )
        ax_conc.plot(times, raw_means, "o-", color="#2c7fb8",
                     markersize=4, linewidth=1, zorder=3, label="Raw mean")
        ax_conc.plot(times, smoothed, color="#e34a33",
                     linewidth=2, zorder=4, label="Smoothed")

        for d_start, d_end in DARK_PERIODS:
            ax_conc.axvspan(d_start, d_end, color="gray", alpha=0.15, zorder=0)
        ax_conc.set_title(met, fontsize=9)
        ax_conc.tick_params(labelsize=8, labelbottom=False)
        ax_conc.set_xlim(-1, 47)
        if col == 0:
            ax_conc.set_ylabel("Conc (nM)", fontsize=9)

        # --- Rate panel ---
        density = np.array([cd_lookup[t] for t in times])
        rate_raw = _rate_from_concs(times, raw_means, density)
        rate_smooth = _rate_from_concs(times, smoothed, density)
        midpoints = (times[:-1] + times[1:]) / 2

        ax_rate.plot(midpoints, rate_raw, "o-", color="#2c7fb8",
                     markersize=3, linewidth=1, label="From raw means")
        ax_rate.plot(midpoints, rate_smooth, "s-", color="#e34a33",
                     markersize=4, linewidth=2, label="From smoothed (used)")
        ax_rate.axhline(0, color="black", linewidth=0.5)

        for d_start, d_end in DARK_PERIODS:
            ax_rate.axvspan(d_start, d_end, color="gray", alpha=0.15, zorder=0)
        ax_rate.tick_params(labelsize=8)
        ax_rate.set_xlim(-1, 47)
        if col == 0:
            ax_rate.set_ylabel("Rate\n(nmol/cell/hr)", fontsize=9)
        if met_row == n_met_rows - 1:
            ax_rate.set_xlabel("Time (h)", fontsize=9)

    # Hide unused axes (slots beyond the number of metabolites)
    for i in range(n_mets, n_met_rows * ncols):
        met_row = i // ncols
        col = i % ncols
        axes[met_row * 2, col].set_visible(False)
        axes[met_row * 2 + 1, col].set_visible(False)

    # Shared legend at the bottom
    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="#666", linestyle="",
                   alpha=0.5, label="Individual replicate measurement"),
        plt.Line2D([0], [0], marker="o", color="#2c7fb8",
                   linestyle="-", label="Raw mean / rate from raw means"),
        plt.Line2D([0], [0], marker="s", color="#e34a33",
                   linestyle="-", linewidth=2,
                   label="Smoothed / rate from smoothed (used downstream)"),
        Patch(facecolor="gray", alpha=0.3, label="Dark period"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        fontsize=9,
        bbox_to_anchor=(0.5, -0.01),
    )

    plt.tight_layout(rect=(0.0, 0.03, 1, 1))
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

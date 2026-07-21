from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

FILE_PATH = Path(__file__).resolve().parent
RES_PATH = FILE_PATH / "results"
FIG_PATH = FILE_PATH / "figures"

# Make the results and figures directories if they don't exist
RES_PATH.mkdir(exist_ok=True)
FIG_PATH.mkdir(exist_ok=True)

# Define colors to use
CUE_C = "#4f8a8b"
GROWTH_C = "#e8896a"


def main():
    # Load the results
    results_df = pd.read_csv(RES_PATH / "growth_and_cue.csv")
    substrate_df = pd.read_csv(RES_PATH / "substrate_panel.csv")
    # Rename the "name" column in the substrate_df to "substrate" so that it can be merged with the results_df
    substrate_df = substrate_df.rename(columns={"name": "substrate"})
    # Merge the results with the substrate info
    merged_df = pd.merge(results_df, substrate_df, on="met_id")

    # For debugging
    print(merged_df.groupby("entry_point")[["atp_cost", "growth_rate"]].mean())

    # Fit a regression line for growth rate and atp_cost
    m_growth_v_atp, d_growth_v_atp = fit(merged_df, "growth_rate", ["atp_cost"])
    m_cue_v_atp, d_cue_v_atp = fit(merged_df, "cue", ["atp_cost"])

    # Fit a regression line for CUE and NOSC
    m_cue, d_cue = fit(merged_df, "cue", ["nosc"])
    # Fit a regression line for growth rate and NOSC
    m_growth, d_growth = fit(merged_df, "growth_rate", ["nosc"])

    # Plot the results on twin axes
    fig, ax = plt.subplots(figsize=(5, 5))
    ax2 = ax.twinx()
    reg_plot(ax, m_cue, d_cue, "cue", ["nosc"], colors=CUE_C)
    reg_plot(ax2, m_growth, d_growth, "growth_rate", ["nosc"], colors=GROWTH_C)
    # Set the Y axis labels and include the R2 values in the label and color to match the points
    ax.set_ylabel(
        f"CUE (R² = {m_cue.score(d_cue[['nosc']], d_cue['cue']):.2f})", color=CUE_C
    )
    ax2.set_ylabel(
        f"growth rate (R² = {m_growth.score(d_growth[['nosc']], d_growth['growth_rate']):.2f})",
        color=GROWTH_C,
    )
    # Color the tick lables to match the points/axis
    ax.tick_params(axis="y", colors=CUE_C)
    ax2.tick_params(axis="y", colors=GROWTH_C)

    plt.tight_layout()
    plt.savefig(FIG_PATH / "nosc_regression.png")
    # Save as an editedable format (e.g. SVG)
    plt.savefig(FIG_PATH / "nosc_regression.svg")

    # Make a table of the R2 values
    r2_df = pd.DataFrame(
        {
            "metric": ["CUE", "growth rate"],
            "NOSC": [
                m_cue.score(d_cue[["nosc"]], d_cue["cue"]),
                m_growth.score(d_growth[["nosc"]], d_growth["growth_rate"]),
            ],
            "ATP cost": [
                m_cue_v_atp.score(d_cue_v_atp[["atp_cost"]], d_cue_v_atp["cue"]),
                m_growth_v_atp.score(
                    d_growth_v_atp[["atp_cost"]], d_growth_v_atp["growth_rate"]
                ),
            ],
        }
    )
    # Round everything to 3 decimal places
    r2_df["NOSC"] = r2_df["NOSC"].round(3)
    r2_df["ATP cost"] = r2_df["ATP cost"].round(3)
    r2_df.to_csv(RES_PATH / "regression_r2s.csv", index=False)


# Helper function to fit a regression model and return the model and the dataframe used for fitting
def fit(df, y, xcols):
    d = df.dropna(subset=[y] + xcols)
    m = LinearRegression().fit(d[xcols], d[y])
    return m, d


# Helper function to make a regression plot
def reg_plot(ax, m, d, y, xcols, colors=None):
    if len(xcols) == 1:  # single predictor: observed vs the predictor
        x = d[xcols[0]]
        xlab = xcols[0]
    else:  # multiple: observed vs fitted
        x = m.predict(d[xcols].values)
        xlab = "predicted " + y

    ax.scatter(x, d[y], c=colors, s=70, edgecolor="white", linewidth=0.8, zorder=3)
    xs = sorted(x)  # regression/1:1 line
    if len(xcols) == 1:
        ax.plot(
            xs,
            m.intercept_ + m.coef_[0] * pd.Series(xs),
            "--",
            color="#666",
            lw=1,
        )
    else:
        ax.plot(xs, xs, "--", color="#666", lw=1)
    ax.set_xlabel(xlab)
    ax.set_ylabel(y)
    ax.spines["top"].set_visible(False)


if __name__ == "__main__":
    main()

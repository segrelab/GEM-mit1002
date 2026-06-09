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


def main():
    # Load the results
    results_df = pd.read_csv(RES_PATH / "growth_and_cue.csv")
    substrate_df = pd.read_csv(RES_PATH / "substrate_panel.csv")
    # Rename the "name" column in the substrate_df to "substrate" so that it can be merged with the results_df
    substrate_df = substrate_df.rename(columns={"name": "substrate"})
    # Merge the results with the substrate info
    merged_df = pd.merge(results_df, substrate_df, on="met_id")

    # Fit a regression line for CUE and NOSC
    m_cue, d_cue = fit(merged_df, "cue", ["nosc"])

    # Print the results
    print(f"CUE ~ nosc: R² = {m_cue.score(d_cue[['nosc']], d_cue['cue']):.3f}")

    # Plot the results
    fig, ax = plt.subplots(figsize=(5, 5))
    reg_plot(ax, m_cue, d_cue, "cue", ["nosc"])
    plt.tight_layout()
    plt.savefig(FIG_PATH / "cue_vs_nosc.png")


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
    ax.text(
        0.05,
        0.9,
        f"$R^2$={m.score(d[xcols], d[y]):.2f}",
        transform=ax.transAxes,
        fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)


if __name__ == "__main__":
    main()

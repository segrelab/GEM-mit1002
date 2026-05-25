"""Load the results from the PCA generated in run_pca.py and plot:
- 2D scatter of the the conditions in PC space
- Bar plot of the top contributing reactions to each PC"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)


def main():
    # Load the PCA scores and loadings
    scores_df = pd.read_csv(OUT_PATH / "pca_scores.csv", index_col=0)
    loadings_df = pd.read_csv(OUT_PATH / "pca_loadings.csv", index_col=0)

    # Load the explained variance ratios from the PCA run
    explained_var = pd.read_csv(OUT_PATH / "explained_variance.csv", index_col=0)[
        "explained_variance_ratio"
    ].values

    # Plot the scores (conditions in PC space)
    plot_scores(scores_df, explained_var)

    # Plot the top contributing reactions to each PC
    plot_top_loadings(
        loadings_df, top_n=15, model=None
    )  # Pass model if you want nicer labels


def plot_scores(scores_df, explained_var):
    fig, ax = plt.subplots(figsize=(8, 7))

    # Color points by substrate. Use the same palette you've been using elsewhere
    # so the colors are consistent across figures.
    colors = sns.color_palette("tab10", n_colors=len(scores_df))

    for (substrate, row), color in zip(scores_df.iterrows(), colors):
        ax.scatter(
            row["PC1"],
            row["PC2"],
            color=color,
            s=200,
            edgecolor="black",
            linewidth=0.5,
        )
        ax.annotate(
            substrate,
            (row["PC1"], row["PC2"]),
            xytext=(8, 0),
            textcoords="offset points",
            fontsize=10,
            color="gray",
            va="center",
        )

    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.set_xlabel(f"PC1 ({explained_var[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({explained_var[1]:.1%} variance)")
    ax.set_title("PCA of MIT1002 flux distributions across substrates")
    plt.tight_layout()
    fig.savefig(OUT_PATH / "pca_scores.png", dpi=150, bbox_inches="tight")


def plot_top_loadings(loadings_df, top_n=15, model=None):
    n_pcs_to_show = 2  # PC1 and PC2
    fig, axes = plt.subplots(1, n_pcs_to_show, figsize=(14, max(6, 0.35 * top_n)))

    for i, pc in enumerate([f"PC{j+1}" for j in range(n_pcs_to_show)]):
        # Get top N reactions by absolute loading
        top_reactions = loadings_df[pc].abs().nlargest(top_n).index
        # Get the actual (signed) loadings
        loadings_to_plot = loadings_df.loc[top_reactions, pc].sort_values()

        # Build informative labels if model is available
        if model is not None:
            labels = []
            for rxn_id in loadings_to_plot.index:
                try:
                    rxn = model.reactions.get_by_id(rxn_id)
                    name = rxn.name[:30] if rxn.name else rxn_id
                    labels.append(f"{rxn_id}: {name}")
                except KeyError:
                    labels.append(rxn_id)
        else:
            labels = loadings_to_plot.index.tolist()

        colors = ["#d95f0e" if x < 0 else "#2c7fb8" for x in loadings_to_plot.values]
        axes[i].barh(
            range(len(loadings_to_plot)), loadings_to_plot.values, color=colors
        )
        axes[i].set_yticks(range(len(loadings_to_plot)))
        axes[i].set_yticklabels(labels, fontsize=8)
        axes[i].axvline(0, color="black", linewidth=0.5)
        axes[i].set_xlabel(f"{pc} loading")
        axes[i].set_title(f"Top {top_n} reactions by absolute loading on {pc}")

    plt.tight_layout()
    plt.savefig(OUT_PATH / "pca_loadings_bars.png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()

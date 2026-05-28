"""Make the 9x9 pairwise epistasis heatmap."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Load the original metabolite list to get their concentrations to sort by
TOP_10_PATH = TOP_10_DIR / "data" / "top10_exometabolites.csv"
top_10_exometabolites = pd.read_csv(TOP_10_PATH)
# Filter the data to only include the Prochlorococcus marinus metabolites
top_10_exometabolites = top_10_exometabolites[
    top_10_exometabolites["organism"] == "Prochlorococcus marinus"
].copy()
# Remove the metabolites that don't support growth/that I didn't test
# In other scripts I use the model to check this, but I don't want to have to
# load the model here just for that, and I know from before that it's just Met
# and Phe that don't have exchanges
top_10_exometabolites = top_10_exometabolites[
    ~top_10_exometabolites["metabolite"].isin(["methionine", "phenylalanine"])
].copy()
# Sort that list by the carbon concentration
top_10_exometabolites = top_10_exometabolites.sort_values(
    by="carbon_concentration", ascending=False
).copy()

# Get the desired order of substrates
substrate_order = top_10_exometabolites["metabolite"].tolist()

# Load the results from the pairwise simulations and the single substrate
# simulations from the single + cocktail simulations
pairs = pd.read_csv(OUT_PATH / "pairwise.csv")
singles = pd.read_csv(
    TOP_10_DIR
    / "single_and_cocktail_sims"
    / "results"
    / "single_and_cocktail_results.csv"
)
singles = singles[singles["condition"] == "single"]
singles_lookup = dict(zip(singles["substrate"], singles["growth_rate"]))

# For each pair, compute the epistasis score:
# epistasis = growth(A+B) - (growth(A) + growth(B)) / 2
# (The /2 is because in the pair, each substrate is at half its single dose;
#  alternative is to compare to max or to a specific null model)
pairs["expected"] = pairs.apply(
    lambda r: (singles_lookup[r["substrate_a"]] + singles_lookup[r["substrate_b"]]) / 2,
    axis=1,
)
pairs["epistasis"] = pairs["growth_rate"] - pairs["expected"]

# Reshape into 9x9 matrix
# TODO: Since I reorder based on the clustering, I don't need to use the
# substrate order here, that would help me remove unneccesarry code loading and
# sorting the metabolite data
matrix = pd.DataFrame(index=substrate_order, columns=substrate_order, dtype=float)
for _, row in pairs.iterrows():
    matrix.loc[row["substrate_a"], row["substrate_b"]] = row["epistasis"]
    matrix.loc[row["substrate_b"], row["substrate_a"]] = row["epistasis"]

# Round anything with an absolute value less than 0.001 to 0, to make the plot cleaner
matrix = matrix.map(lambda x: 0 if abs(x) < 0.001 else x)

# Fill the diagnonal with 0s (the epistasis of a substrate with itself is 0)
np.fill_diagonal(matrix.values, 0)

# Save the matrix to a csv file
matrix.to_csv(OUT_PATH / "epistasis_matrix.csv")

# Create a custom sequential colormap
# TODO: Acuatlly import plot_styles, rather than copying a hexcode
custom_cmap = sns.light_palette("#024064", as_cmap=True)

# Plot with colormap that begins at 0
# (positive = synergy, negative = antagonism)
g = sns.clustermap(
    matrix,
    cmap=custom_cmap,
    vmin=0,
    annot=True,  # Show the values in the cells
    fmt=".3f",  # Format annotations to 3 decimal places
    linewidths=0.5,
    method="average",  # Clustering method
    metric="euclidean",  # Distance metric for clustering
    figsize=(10, 10),
    cbar_kws={"label": "Epistasis Score (Growth Rate Difference)"},
)

# Style the axes
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", color="gray"
)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), color="gray")

# Move the colorbar
# Get the heatmap's position
heatmap_pos = g.ax_heatmap.get_position()
# Move the colorbar to be flush with the right side of the heatmap, full height
g.cax.set_position(
    [
        heatmap_pos.x1 + 0.12,  # left: 10% right of heatmap's right edge
        heatmap_pos.y0,  # bottom: aligned with heatmap's bottom
        0.02,  # width: 2% of figure
        heatmap_pos.height,  # height: matches heatmap exactly
    ]
)

# Change the colorbar title, ticks, and tick labels to be gray
g.cax.tick_params(colors="gray")  # Set the ticks and tick labels
g.cax.yaxis.label.set_color("gray")  # Set the label text

g.savefig(OUT_PATH / "figure_epistasis.png", dpi=150, bbox_inches="tight")

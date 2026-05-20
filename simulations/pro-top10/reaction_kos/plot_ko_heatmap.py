"""Run all single reaction knockouts on all substrates (single and cocktail)
and save the results."""

from pathlib import Path

import cobra
import numpy as np
import pickle as pkl
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Set a minimum growth rate threshold for essentiality
MIN_GROWTH_RATE = 0.0001  # 1 / hr

# Define the maintenance reaction ID
MAINTENANCE_REACTION_ID = "rxn00062_c0"

# Load the model
MODEL = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

# Load the ko results
ko_results_df = pd.read_csv(OUT_PATH / "single_reaction_ko_results.csv", index_col=0)

# Round the growth rates below the minimum growth rate threshold to 0 to make the heatmap easier to read
ko_results_df = ko_results_df.map(lambda x: 0 if x < MIN_GROWTH_RATE else x)

# Filter the data to only include reactions that are essential in at least one
# condition (growth rate < minimum growth rate threshold), but not all
semi_essential_reactions_df = ko_results_df[
    ko_results_df.lt(MIN_GROWTH_RATE).any(axis=1)
    & ~ko_results_df.lt(MIN_GROWTH_RATE).all(axis=1)
]

# Normalize the growth rates by the growth rate of the wild type on each substrate
# To get the wt growth rate, find the max of the column (excluding for the maintenance reaction KO)
# This is because the maintenance reaction KO will have a higher growth rate than the wt
wt_growth_rates = ko_results_df.drop(MAINTENANCE_REACTION_ID, axis=0).max()
semi_essential_reactions_df_norm = semi_essential_reactions_df.div(wt_growth_rates)

# Plot a heatmap of the growth rates for each reaction knockout and substrate
plt.figure(figsize=(10, 20))
sns.heatmap(
    semi_essential_reactions_df_norm,
    cmap="viridis",
    cbar_kws={"label": "Growth rate (1/hr)"},
)
plt.xlabel("Substrate")
plt.ylabel("Knocked out reaction")
plt.title(
    "Growth rates for single reaction knockouts on single substrates and the cocktail"
)
plt.tight_layout()

# Save the plot
plt.savefig(OUT_PATH / "ko_heatmap.png", dpi=300)

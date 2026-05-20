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

# Pathway groupings
# TODO: Double check these groupings and make sure they are correct
PATHWAY_GROUPS = {
    "rxn15021_c0": "Val biosynthesis",
    "rxn15466_c0": "Val biosynthesis",
    "rxn15467_c0": "Val biosynthesis",
    "rxn03435_c0": "Ile biosynthesis",
    "rxn03436_c0": "Ile biosynthesis",
    "rxn03437_c0": "Ile biosynthesis",
    "rxn08043_c0": "Ile biosynthesis",
    "rxn00737_c0": "Ile biosynthesis",
    "rxn00902_c0": "Leu biosynthesis",
    "rxn01208_c0": "Leu biosynthesis",
    "rxn02789_c0": "Leu biosynthesis",
    "rxn02811_c0": "Leu biosynthesis",
    "rxn03062_c0": "Leu biosynthesis",
    "rxn00974_c0": "Glyoxylate / TCA",
    "rxn01388_c0": "Glyoxylate / TCA",
    "rxn00256_c0": "Glyoxylate / TCA",
    "rxn00336_c0": "Glyoxylate / TCA",
    "rxn09272_c0": "TCA (SDH)",
    "rxn00692_c0": "Glycine cleavage / C1",
    "rxn09498_c0": "Glycine cleavage / C1",
    "rxn09499_c0": "Glycine cleavage / C1",
    "rxn01241_c0": "Glycine cleavage / C1",
    "rxn05468_c0": "O2 transport",
    "EX_cpd00007_e0": "O2 transport",
    "rxn00102_c0": "Other",
}

# Define a deterministic pathway display order
PATHWAY_ORDER = [
    "Val biosynthesis",
    "Ile biosynthesis",
    "Leu biosynthesis",
    "Glyoxylate / TCA",
    "TCA (SDH)",
    "Glycine cleavage / C1",
    "O2 transport",
    "Other",
]

# Load the model
MODEL = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")


def build_label(rxn_id, model):
    """Build a y-axis label like 'rxn00974_c0: aconitase: Citrate <=> ...'."""
    try:
        rxn = model.reactions.get_by_id(rxn_id)
        # Truncate long reaction strings so labels don't dominate the figure
        rxn_str = rxn.build_reaction_string(use_metabolite_names=True)
        if len(rxn_str) > 45:
            rxn_str = rxn_str[:42] + "..."
        name = rxn.name if rxn.name else ""
        if len(name) > 30:
            name = name[:27] + "..."
        return f"{rxn_id}: {name} | {rxn_str}"
    except KeyError:
        return rxn_id


# Load the ko results
ko_results_df = pd.read_csv(OUT_PATH / "single_reaction_ko_results.csv", index_col=0)

# Round the growth rates below the minimum growth rate threshold to 0 to make the heatmap easier to read
ko_results_df = ko_results_df.map(lambda x: 0 if x < MIN_GROWTH_RATE else x)

# Filter the data to only include reactions that are essential in at least one
# condition (growth rate < minimum growth rate threshold), but not all
cond_df = ko_results_df[
    ko_results_df.lt(MIN_GROWTH_RATE).any(axis=1)
    & ~ko_results_df.lt(MIN_GROWTH_RATE).all(axis=1)
]

# Normalize the growth rates by the growth rate of the wild type on each substrate
# To get the wt growth rate, find the max of the column (excluding for the maintenance reaction KO)
# This is because the maintenance reaction KO will have a higher growth rate than the wt
wt_growth_rates = ko_results_df.drop(MAINTENANCE_REACTION_ID, axis=0).max()
cond_df_norm = cond_df.div(wt_growth_rates)

# Sort rows by pathway, then by reaction ID within pathway
cond_df_norm = cond_df_norm.copy()
cond_df_norm["pathway"] = cond_df_norm.index.map(
    lambda x: PATHWAY_GROUPS.get(x, "Unassigned")
)
cond_df_norm["pathway_order"] = cond_df_norm["pathway"].map(
    {p: i for i, p in enumerate(PATHWAY_ORDER)}
)
# Sort and then drop the helper columns
cond_df_norm = cond_df_norm.sort_values(["pathway_order", "pathway"])
pathways = cond_df_norm["pathway"].tolist()
cond_df_norm = cond_df_norm.drop(columns=["pathway", "pathway_order"])

# Build informative labels
labels = [build_label(rxn_id, MODEL) for rxn_id in cond_df_norm.index]

# Plot
n_rows = len(cond_df_norm)
fig = plt.figure(figsize=(10, max(6, 0.4 * n_rows)))

# Define axes positions manually using add_axes
# [left, bottom, width, height] in figure coordinates
heatmap_ax = fig.add_axes([0.35, 0.1, 0.4, 0.85])  # heatmap occupies center
cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.5])  # colorbar far right

sns.heatmap(
    cond_df_norm.values,
    cmap="viridis",
    cbar_kws={"label": "Growth rate (normalized to WT)"},
    yticklabels=labels,
    xticklabels=cond_df_norm.columns,
    linewidths=0.3,
    linecolor="lightgray",
    ax=heatmap_ax,
    cbar_ax=cbar_ax,
)

# Add pathway group separators
for i in range(1, len(pathways)):
    if pathways[i] != pathways[i - 1]:
        heatmap_ax.axhline(i, color="red", linewidth=2)

# Add pathway labels on the right side
seen = set()
for i, p in enumerate(pathways):
    if p in seen:
        continue
    # Find the midpoint of this pathway's rows
    idx_in_group = [j for j, x in enumerate(pathways) if x == p]
    mid = (min(idx_in_group) + max(idx_in_group)) / 2 + 0.5
    heatmap_ax.text(
        len(cond_df_norm.columns) + 0.3,
        mid,
        p,
        fontsize=10,
        va="center",
        color="darkred",
    )
    seen.add(p)

heatmap_ax.set_xlabel("Substrate condition")
heatmap_ax.set_ylabel("")
heatmap_ax.set_title("Conditionally essential reactions across Pro top-10 substrate conditions")
plt.xticks(rotation=45, ha="right")

fig.savefig(OUT_PATH / "ko_heatmap.png", dpi=300, bbox_inches="tight")

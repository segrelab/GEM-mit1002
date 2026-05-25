from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"

# Add the repo root to the system path so we can import from the plot_styles file
sys.path.append(str(REPO_ROOT))
import plot_styles  # Import the plot styles from the repo

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Load the original metabolite list to get their concentrations
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

# Load the results from the single substrate and cocktail simulations
results = pd.read_csv(
    TOP_10_DIR
    / "single_and_cocktail_sims"
    / "results"
    / "single_and_cocktail_results.csv"
)

# Calculate the sum of parts using pandas operations
# Filter for single substrate results
singles_df = results[results["condition"] == "single"].copy()
# Calculate the weight for each metabolite
top_10_exometabolites["weight"] = (
    top_10_exometabolites["carbon_concentration"]
    / top_10_exometabolites["carbon_concentration"].sum()
)
# Merge the growth rates with their weights
merged_df = pd.merge(
    singles_df,
    top_10_exometabolites[["metabolite", "weight"]],
    left_on="substrate",
    right_on="metabolite",
)
# Calculate the weighted growth rate for each single substrate
merged_df["weighted_growth"] = merged_df["growth_rate"] * merged_df["weight"]
# Compute the final sum of parts
sum_of_parts = merged_df["weighted_growth"].sum()
# Add the sum of parts to the results DataFrame
sum_row = pd.DataFrame(
    [
        {
            "condition": "sum_of_parts",
            "substrate": "sum_of_parts",
            "growth_rate": sum_of_parts,
        }
    ]
)
results = pd.concat([results, sum_row], ignore_index=True)

# Select the colors I want to use from the color palette
colors = [
    plot_styles.tsitp_colors["dark_blue"],
    plot_styles.tsitp_colors["dark_green"],
    plot_styles.tsitp_colors["dark_orange"],
]

# Plot the results
fig, ax = plt.subplots(figsize=(8, 6))
sns.barplot(
    x="substrate",
    y="growth_rate",
    data=results,
    hue="condition",
    ax=ax,
    palette=colors,
)
plt.xlabel("Substrate")
plt.ylabel("Growth Rate (1/hr)")
plt.title("Growth Rate on Single Substrates and Cocktail")
plt.xticks(rotation=45)

# Move the legend to the right of the plot
ax.legend(loc="upper center", bbox_to_anchor=(1.15, 0.5), ncol=1, frameon=False)

plt.tight_layout()

# Apply the plot style
plot_styles.set_plot_style(ax)

plt.savefig(OUT_PATH / "growth_rate_comparison.png")

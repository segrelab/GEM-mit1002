"""Make a stacked barchart of the exudation fluxes on the single substrates and
the cocktail."""

from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Load the model
# So I can get metabolite names from reaction IDs
model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

# Load the results from the simulations from the single + cocktail simulations
# Read the "fluxes" column as a dictionary
results = pd.read_csv(
    TOP_10_DIR
    / "single_and_cocktail_sims"
    / "results"
    / "single_and_cocktail_results.csv",
    converters={"fluxes": eval},
)

# Extract the exchange fluxes from the 'fluxes' column
ex_fluxes_list = results["fluxes"].apply(
    lambda d: {k: v for k, v in d.items() if k.startswith("EX_")}
)

# Create a new DataFrame from the list of dictionaries
ex_fluxes_df = pd.DataFrame(ex_fluxes_list.tolist(), index=results["substrate"])

# There are many columns (reactions) that have all 0s
# Drop those so that the legend is smaller/the plot uses fewer colors
# Drop all columns that have all 0s
ex_fluxes_df = ex_fluxes_df.loc[:, (ex_fluxes_df != 0).any(axis=0)]

# Rename the columns with the name of the metabolite instead of the reaction ID
new_column_names = {}
for rxn_id in ex_fluxes_df.columns:
    rxn = model.reactions.get_by_id(rxn_id)
    # Get the corresponding metabolite from the reaction
    # The exchange reactions have only one metabolite, so get that one
    met = list(rxn.metabolites.keys())[0]
    # Get the name of that metabolite
    met_name = met.name
    new_column_names[rxn_id] = met_name
# Rename the columns with the new names
ex_fluxes_df = ex_fluxes_df.rename(columns=new_column_names)

# Save the dataframe to a CSV file
ex_fluxes_df.to_csv(OUT_PATH / "exchange_fluxes_single_and_cocktail.csv")

# Make a copy of the dataframe with only the positive fluxes (only exudation)
exudation_df = ex_fluxes_df.map(lambda x: x if x > 0 else 0)
# Drop all columns that have all 0s
exudation_df = exudation_df.loc[:, (exudation_df != 0).any(axis=0)]

# Make a copy of the dataframe with only the negative fluxes (only uptake)
uptake_df = ex_fluxes_df.map(lambda x: x if x < 0 else 0)
# Drop all columns that have all 0s
uptake_df = uptake_df.loc[:, (uptake_df != 0).any(axis=0)]

# Plot both the exudation and uptake fluxes
ex_fluxes_df.plot(kind="bar", stacked=True, figsize=(12, 7))
# Add labels and title for clarity
plt.title("Exchange Fluxes")
plt.xlabel("Solution")
plt.ylabel("Reaction Flux (mmol/gDW/hr)")
plt.legend(title="Exchange Reactions", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
# Save the plot
plt.savefig(OUT_PATH / "exchange_fluxes_single_and_cocktail.png")

# Plot the exudation fluxes
exudation_df.plot(kind="bar", stacked=True, figsize=(12, 7))
# Add labels and title for clarity
plt.title("Exudation Fluxes")
plt.xlabel("Solution")
plt.ylabel("Reaction Flux (mmol/gDW/hr)")
plt.legend(title="Metabolites Released", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
# Save the plot
plt.savefig(OUT_PATH / "exudation_fluxes_single_and_cocktail.png")

# Plot the uptake fluxes
uptake_df.plot(kind="bar", stacked=True, figsize=(12, 7))
# Add labels and title for clarity
plt.title("Uptake Fluxes")
plt.xlabel("Solution")
plt.ylabel("Reaction Flux (mmol/gDW/hr)")
plt.legend(title="Metabolites Taken Up", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
# Save the plot
plt.savefig(OUT_PATH / "uptake_fluxes_single_and_cocktail.png")

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


# Helper function to apply to each row
def determine_direction(column, e_met):
    # Get the reaction ID (the column name)
    rxn_id = column.name
    # Get the reaction object
    rxn = model.reactions.get_by_id(rxn_id)
    # Check if the extracellular metabolite is a reactant or product in the reaction
    if e_met in rxn.reactants:
        # If the flux is positive, the metabolite is moving into the cell, so return "in"
        if flux > 0:
            return "in"
        # If the flux is negative, the metabolite is moving out of the cell, so return "out"
        elif flux < 0:
            return "out"
    elif e_met in rxn.products:
        # If the flux is positive, the metabolite is moving out of the cell, so return "out"
        if flux > 0:
            return "out"
        # If the flux is negative, the metabolite is moving into the cell, so return "in"
        elif flux < 0:
            return "in"


# Load the model
# So I can get metabolite names from reaction IDs
model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

# Get the extracellular proton and sodium metabolites
h_e_met = model.metabolites.cpd00067_e0
na_e_met = model.metabolites.cpd00971_e0

# Load the results from the simulations from the single + cocktail simulations
# Read the "fluxes" column as a dictionary
results = pd.read_csv(
    TOP_10_DIR
    / "single_and_cocktail_sims"
    / "results"
    / "single_and_cocktail_results.csv",
    converters={"fluxes": eval},
)

# Get a list of reaction IDs that involve extracellular protons
h_e_rxns = [r.id for r in model.reactions if h_e_met in r.metabolites]
# Get a list of reaction IDs that involve extracellular sodium
na_e_rxns = [r.id for r in model.reactions if na_e_met in r.metabolites]

# Extract the exchange fluxes from the 'fluxes' column
na_fluxes_list = results["fluxes"].apply(
    lambda d: {k: v for k, v in d.items() if k in na_e_rxns}
)

# Create a new DataFrame from the list of dictionaries
na_fluxes_df = pd.DataFrame(na_fluxes_list.tolist(), index=results["substrate"])

# There are many columns (reactions) that have all 0s
# Drop those so that the legend is smaller/the plot uses fewer colors
# Drop all columns that have all 0s
na_fluxes_df = na_fluxes_df.loc[:, (na_fluxes_df != 0).any(axis=0)]

# Make the sign match the direction of Na movement
# Positive means Na is moving out of the cell, negative means Na is moving into the cell
# This also depends on the "side" of the reaction that the extracellular sodium is on
# If the extracellular sodium is a reactant and the flux is positive, the sodium is moving into the cell, so make it negative
# If the extracellular sodium is a reactant and the flux is negative, the sodium is moving out of the cell, so make it positive
# If the extracellular sodium is a product and the flux is positive, the sodium is moving out of the cell, so keep it positive
# If the extracellular sodium is a product and the flux is negative, the sodium is moving into the cell, so keep it negative
# We can simplify all of this to if extracellular sodium is a reactant, multiply the flux by -1, if it's a product, keep it the same
for rxn_id in na_fluxes_df.columns:
    rxn = model.reactions.get_by_id(rxn_id)
    if na_e_met in rxn.reactants:
        na_fluxes_df[rxn_id] = na_fluxes_df[rxn_id].apply(lambda flux: -flux)

# Rename the columns with the name of the reaction string of the reaction ID
new_column_names = {}
for rxn_id in na_fluxes_df.columns:
    # Get the reaction
    rxn = model.reactions.get_by_id(rxn_id)
    # Build the reaction string
    reaction_string = rxn.build_reaction_string(use_metabolite_names=True)
    # Set the new column name to be the reaction string
    new_column_names[rxn_id] = reaction_string
# Rename the columns with the new names
na_fluxes_df = na_fluxes_df.rename(columns=new_column_names)

# Save the dataframe to a CSV file
na_fluxes_df.to_csv(OUT_PATH / "smf_fluxes_single_and_cocktail.csv")

# # Make a copy of the dataframe with only the positive fluxes (only exudation)
# exudation_df = na_fluxes_df.map(lambda x: x if x > 0 else 0)
# # Drop all columns that have all 0s
# exudation_df = exudation_df.loc[:, (exudation_df != 0).any(axis=0)]

# # Make a copy of the dataframe with only the negative fluxes (only uptake)
# uptake_df = na_fluxes_df.map(lambda x: x if x < 0 else 0)
# # Drop all columns that have all 0s
# uptake_df = uptake_df.loc[:, (uptake_df != 0).any(axis=0)]

# Plot both the exudation and uptake fluxes
na_fluxes_df.plot(kind="bar", stacked=True, figsize=(12, 7))
# Add labels and title for clarity
plt.title("Exchange Fluxes")
plt.xlabel("Solution")
plt.ylabel("Reaction Flux (mmol/gDW/hr)")
plt.legend(title="Exchange Reactions", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
# Save the plot
plt.savefig(OUT_PATH / "exchange_fluxes_single_and_cocktail.png")

# # Plot the exudation fluxes
# exudation_df.plot(kind="bar", stacked=True, figsize=(12, 7))
# # Add labels and title for clarity
# plt.title("Exudation Fluxes")
# plt.xlabel("Solution")
# plt.ylabel("Reaction Flux (mmol/gDW/hr)")
# plt.legend(title="Metabolites Released", bbox_to_anchor=(1.05, 1), loc="upper left")
# plt.tight_layout()
# # Save the plot
# plt.savefig(OUT_PATH / "exudation_fluxes_single_and_cocktail.png")

# # Plot the uptake fluxes
# uptake_df.plot(kind="bar", stacked=True, figsize=(12, 7))
# # Add labels and title for clarity
# plt.title("Uptake Fluxes")
# plt.xlabel("Solution")
# plt.ylabel("Reaction Flux (mmol/gDW/hr)")
# plt.legend(title="Metabolites Taken Up", bbox_to_anchor=(1.05, 1), loc="upper left")
# plt.tight_layout()
# # Save the plot
# plt.savefig(OUT_PATH / "uptake_fluxes_single_and_cocktail.png")

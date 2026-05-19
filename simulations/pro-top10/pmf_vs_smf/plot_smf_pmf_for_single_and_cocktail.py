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
MODEL = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

# Load the results from the simulations from the single + cocktail simulations
# Read the "fluxes" column as a dictionary
RESULTS = pd.read_csv(
    TOP_10_DIR
    / "single_and_cocktail_sims"
    / "results"
    / "single_and_cocktail_results.csv",
    converters={"fluxes": eval},
)


def main():
    # Get the extracellular proton and sodium metabolites
    h_e_met = MODEL.metabolites.cpd00067_e0
    na_e_met = MODEL.metabolites.cpd00971_e0

    # Do for sodium
    plot_fluxes_for_pumped_met(na_e_met, "smf")
    # Do for protins
    plot_fluxes_for_pumped_met(h_e_met, "pmf")


def plot_fluxes_for_pumped_met(met, file_prefix):
    # Get a list of reaction IDs that involve the specified metabolite
    met_rxns = [r.id for r in MODEL.reactions if met in r.metabolites]

    # Extract the relevant fluxes from the 'fluxes' column
    met_fluxes_list = RESULTS["fluxes"].apply(
        lambda d: {k: v for k, v in d.items() if k in met_rxns}
    )

    # Create a new DataFrame from the list of dictionaries
    met_fluxes_df = pd.DataFrame(met_fluxes_list.tolist(), index=RESULTS["substrate"])

    # There are many columns (reactions) that have all 0s
    # Drop those so that the legend is smaller/the plot uses fewer colors
    # Drop all columns that have all 0s
    met_fluxes_df = met_fluxes_df.loc[:, (met_fluxes_df != 0).any(axis=0)]

    # Make the sign match the direction of the specified metabolite movement
    # Positive means the metabolite is moving out of the cell, negative means it is moving into the cell
    # This also depends on the "side" of the reaction that the extracellular metabolite is on
    # If the extracellular metabolite is a reactant and the flux is positive, it is moving into the cell, so make it negative
    # If the extracellular metabolite is a reactant and the flux is negative, the metabolite is moving out of the cell, so make it positive
    # If the extracellular metabolite is a product and the flux is positive, the metabolite is moving out of the cell, so keep it positive
    # If the extracellular metabolite is a product and the flux is negative, the metabolite is moving into the cell, so keep it negative
    # We can simplify all of this to if extracellular metabolite is a reactant, multiply the flux by -1, if it's a product, keep it the same
    for rxn_id in met_fluxes_df.columns:
        # Get the reaction object
        rxn = MODEL.reactions.get_by_id(rxn_id)
        # Get the coefficient of the extracellular metabolite in the reaction
        s_coeff = rxn.metabolites[met]
        # Multiply the flux by the coefficient to get the correct rate of metabolite movement
        # Only necessary if the coefficient is not -1 or 1
        if abs(rxn.metabolites[met]) != 1:
            met_fluxes_df[rxn_id] = met_fluxes_df[rxn_id].apply(
                lambda flux: flux * abs(s_coeff)
            )
        # If the extracellular metabolite is a reactant, multiply the flux by -1 to get the correct direction
        if met in rxn.reactants:
            met_fluxes_df[rxn_id] = met_fluxes_df[rxn_id].apply(lambda flux: -flux)

    # Rename the columns with the name of the reaction string of the reaction ID
    new_column_names = {}
    for rxn_id in met_fluxes_df.columns:
        # Get the reaction
        rxn = MODEL.reactions.get_by_id(rxn_id)
        # Build the reaction string
        reaction_string = rxn.build_reaction_string(use_metabolite_names=True)
        # Set the new column name to be the reaction string
        new_column_names[rxn_id] = reaction_string
    # Rename the columns with the new names
    met_fluxes_df = met_fluxes_df.rename(columns=new_column_names)

    # Save the dataframe to a CSV file
    met_fluxes_df.to_csv(OUT_PATH / (file_prefix + "_fluxes_single_and_cocktail.csv"))

    # Plot both the exudation and uptake fluxes
    met_fluxes_df.plot(kind="bar", stacked=True, figsize=(12, 7))
    # Add labels and title for clarity
    plt.title(f"{met.name} Pumping Fluxes")
    plt.xlabel("Solution")
    plt.ylabel("Reaction Flux (mmol/gDW/hr)")
    plt.legend(title="Reactions", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    # Save the plot
    plt.savefig(OUT_PATH / (file_prefix + "_fluxes_single_and_cocktail.png"))


if __name__ == "__main__":
    main()

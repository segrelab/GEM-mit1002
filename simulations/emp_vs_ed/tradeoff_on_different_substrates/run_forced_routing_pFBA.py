import pickle as pkl
from pathlib import Path

import cobra
import pandas as pd
from gem_utilities import media as media_utils
from yaml import warnings

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Save reaction IDs for key reactions
# Define a reaction to indicate flux through the ED pathway
ed_rxn_id = "rxn01477_c0"
# Define a reaction to indicate flux through the EMP
emp_rxn_id = "rxn00558_c0"
# Define the biomass reaction ID
biomass_rxn_id = "bio1_biomass"

# Set the total carbon uptake to use
TOTAL_UPTAKE = 60  # mmol C / gDW / hr

# Load the model
model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

# Load the media definitions
with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
    media_defs = pkl.load(f)

# Load the same substrate panel that was used for the growth and CUE analysis
substrate_df = pd.read_csv(
    REPO_ROOT / "simulations" / "growth_and_cue" / "results" / "substrate_panel.csv"
)

# Subset the substrate panel to only include substrates whose entry points are
# "Glycolysis" or "ED - KDPG"
substrate_df = substrate_df[
    substrate_df["entry_point"].isin(["Glycolysis", "ED - KDPG"])
]

# Make a dictionary to store the results for each substrate
results = {}
# Loop through all of the substrates to test
for _, row in substrate_df.iterrows():
    name = row["name"]
    # Build a media for the current substrate and set it
    media = media_utils.clean_media(model, media_defs["minimal"])
    substrate_uptake_rate = TOTAL_UPTAKE / row["n_c"]
    media[row["exchange_id"]] = substrate_uptake_rate
    model.medium = media
    # Make a dictionary to store the results
    substrate_results = {}
    # Force flux through the EMP pathway and save the growth rate
    with model:
        # Set the lowerbound of the EMP reaction to the uptake rate
        model.reactions.get_by_id(emp_rxn_id).lower_bound = substrate_uptake_rate
        # Run pFBA
        try:
            solution = cobra.flux_analysis.pfba(model)
            # Save the growth rate
            substrate_results["emp"] = solution.fluxes[biomass_rxn_id]
        except:
            print(f"Error occurred while running pFBA with EMP forcing on {name}")
    # Force flux through the ED pathway and save the growth rate
    with model:
        # Set the lowerbound of the ED reaction to the uptake rate
        model.reactions.get_by_id(ed_rxn_id).lower_bound = substrate_uptake_rate
        # Run pFBA
        try:
            solution = cobra.flux_analysis.pfba(model)
            # Save the growth rate
            substrate_results["ed"] = solution.fluxes[biomass_rxn_id]
        except:
            print(f"Error occurred while running pFBA with ED forcing on {name}")
    # Add the results for this substrate to the main results dictionary
    results[name] = substrate_results

# Convert the results to a pandas DataFrame
results_df = pd.DataFrame.from_dict(results, orient="index")
# Save the results to a CSV file
results_df.to_csv(OUT_PATH / "forced_routing_results.csv", index=True)

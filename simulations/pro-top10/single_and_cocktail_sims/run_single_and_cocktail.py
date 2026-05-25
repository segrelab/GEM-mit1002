"""Run pFBA on each single Pro Top 10 Exometabolite, the cocktail of them all,
and save the full results."""

from pathlib import Path

import cobra
import pickle as pkl
import pandas as pd

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Set a total ammount of carbon to take up
# Same number regardless if growing on a single substrate or a mix
# 60 matches a glucose bound of 10
TOTAL_UPTAKE = 60  # mmol C / gDW / hr

# DEFINE THE BIOMASS REACTION ID
BIOMASS_REACTION_ID = "bio1_biomass"

# Load the model
MODEL_PATH = REPO_ROOT / "model.xml"
model = cobra.io.read_sbml_model(MODEL_PATH)

# Load the list of top 10 exometabolites
TOP_10_PATH = TOP_10_DIR / "data" / "top10_exometabolites.csv"
top_10_exometabolites = pd.read_csv(TOP_10_PATH)

# Sort by the carbon concentration
top_10_exometabolites = top_10_exometabolites.sort_values(
    by="carbon_concentration", ascending=False
).copy()

# Filter the data to only include the Prochlorococcus marinus metabolites
top_10_exometabolites = top_10_exometabolites[
    top_10_exometabolites["organism"] == "Prochlorococcus marinus"
].copy()

# Filter the data to remove the metabolites that don't have an exchange
# reaction, and therefore do not support growth
# This should only remove 2 metabolites, methionine and phenylalanine
top_10_exometabolites["exchange_id"] = "EX_" + top_10_exometabolites["met_id"] + "_e0"
top_10_exometabolites = top_10_exometabolites[
    top_10_exometabolites["exchange_id"].isin([r.id for r in model.reactions])
].copy()

# Load the minimal media definition
with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
    media_definitions = pkl.load(f)
minimal_media = media_definitions["minimal"]

# Make a list to hold the results
results = []

# Test each exometabolite in the top10 df at the full TOTAL_UPTAKE
for index, row in top_10_exometabolites.iterrows():
    # Make a copy of the minimal media and add the current substrate to it
    media = minimal_media.copy()
    # The substrate uptake is the total uptake divided by the number of carbons
    # So that every metabolite has the same amount of carbon
    media[row["exchange_id"]] = TOTAL_UPTAKE / row["n_c"]
    # Set the media
    model.medium = media

    # Run pFBA
    solution = cobra.flux_analysis.pfba(model)
    growth = solution.fluxes[BIOMASS_REACTION_ID]

    # Add the results to the list
    results.append(
        {
            "condition": "single",
            "substrate": row["metabolite"],
            "growth_rate": growth,
            "fluxes": solution.fluxes.to_dict(),
        }
    )

# Create a cocktail condition
# All thes substrates are available in their ratio in the top 10 exometabolites
cocktail_media = minimal_media.copy()
for index, row in top_10_exometabolites.iterrows():
    cocktail_media[row["exchange_id"]] = (
        TOTAL_UPTAKE
        * row["carbon_concentration"]
        / top_10_exometabolites["carbon_concentration"].sum()
    ) / row["n_c"]
model.medium = cocktail_media
solution = cobra.flux_analysis.pfba(model)
growth = solution.fluxes[BIOMASS_REACTION_ID]
results.append(
    {
        "condition": "cocktail",
        "substrate": "cocktail",
        "growth_rate": growth,
        "fluxes": solution.fluxes.to_dict(),
    }
)

# Convert the results to a dataframe and save it
results_df = pd.DataFrame(results)
results_df.to_csv(
    OUT_PATH / "single_and_cocktail_results.csv",
    index=False,
)

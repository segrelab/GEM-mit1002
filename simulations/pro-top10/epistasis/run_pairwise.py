"Run pFBA on each possible combination of the top 10 exometabolites (expcept"

"for Met and Phe, which don't support growth)"
from pathlib import Path

import cobra
from itertools import combinations
import pandas as pd
import pickle as pkl

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

# Filter the data to only include the Prochlorococcus marinus metabolites
top_10_exometabolites = top_10_exometabolites[
    top_10_exometabolites["organism"] == "Prochlorococcus marinus"
].copy()

# Sort by the carbon concentration
top_10_exometabolites = top_10_exometabolites.sort_values(
    by="carbon_concentration", ascending=False
).copy()

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

# Loop through all of the combinations of the remaining Top 10 exometabolite df rows
# Loop through all combinations of DataFrame rows using their indices
for index_a, index_b in combinations(top_10_exometabolites.index, 2):
    # Get the full row for each part of the pair
    row_a = top_10_exometabolites.loc[index_a]
    row_b = top_10_exometabolites.loc[index_b]

    # Make a copy of the minimal media and add the current substrate to it
    media = minimal_media.copy()
    # The substrate uptake is half the total uptake divided by the number of carbons
    # So that every metabolite has the same amount of carbon
    media[row_a["exchange_id"]] = TOTAL_UPTAKE / 2 / row_a["n_c"]
    media[row_b["exchange_id"]] = TOTAL_UPTAKE / 2 / row_b["n_c"]
    # Set the media
    model.medium = media

    # Run pFBA
    solution = cobra.flux_analysis.pfba(model)
    growth = solution.fluxes[BIOMASS_REACTION_ID]

    # Add the results to the list
    results.append(
        {
            "substrate_a": row_a["metabolite"],
            "substrate_b": row_b["metabolite"],
            "growth_rate": growth,
        }
    )

# Convert the results to a DataFrame and save it
results_df = pd.DataFrame(results)
results_df.to_csv(OUT_PATH / "pairwise.csv", index=False)

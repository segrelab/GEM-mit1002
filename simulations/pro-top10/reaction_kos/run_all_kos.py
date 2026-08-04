"""Run all single reaction knockouts on all substrates (single and cocktail)
and save the results."""

from pathlib import Path

import cobra
import numpy as np
import pandas as pd

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]

import sys

sys.path.insert(0, str(REPO_ROOT))

from tools.media import MEDIA  # noqa: E402

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Set a total ammount of carbon to take up
# Same number regardless if growing on a single substrate or a mix
# 60 matches a glucose bound of 10
TOTAL_UPTAKE = 60  # mmol C / gDW / hr

# Load the model
MODEL = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")


def main():
    # Load the list of top 10 exometabolites
    TOP_10_PATH = TOP_10_DIR / "data" / "top10_exometabolites.csv"
    top_10_exometabolites = pd.read_csv(TOP_10_PATH)

    # Filter the data to only include the Prochlorococcus marinus metabolites
    top_10_exometabolites = top_10_exometabolites[
        top_10_exometabolites["organism"] == "Prochlorococcus marinus"
    ].copy()

    # Filter the data to remove the metabolites that don't have an exchange
    # reaction, and therefore do not support growth
    # This should only remove 2 metabolites, methionine and phenylalanine
    top_10_exometabolites["exchange_id"] = (
        "EX_" + top_10_exometabolites["met_id"] + "_e0"
    )
    top_10_exometabolites = top_10_exometabolites[
        top_10_exometabolites["exchange_id"].isin([r.id for r in MODEL.reactions])
    ].copy()

    # Load the minimal media definition
    media_definitions = MEDIA
    minimal_media = media_definitions["minimal"]

    # Make a list to hold the results
    results = []

    # Set the media and do the single reaction knockouts for each substrate condition
    for index, row in top_10_exometabolites.iterrows():
        # Print a message to show progress
        print(f"Running single reaction KOs for {row['metabolite']}...")

        # Make a copy of the minimal media and add the current substrate to it
        media = minimal_media.copy()
        # The substrate uptake is the total uptake divided by the number of carbons
        # So that every metabolite has the same amount of carbon
        media[row["exchange_id"]] = TOTAL_UPTAKE / row["n_c"]
        # Set the media
        MODEL.medium = media

        # Run the reaction knockouts
        ko_results = cobra.flux_analysis.single_reaction_deletion(MODEL)

        # Set growth to NaN where status is not 'optimal'
        ko_results["growth"] = ko_results.apply(
            lambda row: row["growth"] if row["status"] == "optimal" else np.nan, axis=1
        )

        # Rename the "growth" column to be the name of the substrate
        ko_results = ko_results.rename(columns={"growth": row["metabolite"]})

        # Convert the 'ids' column from sets to strings
        ko_results["ids"] = ko_results["ids"].apply(lambda x: list(x)[0])

        # Append the results (minus the "status" column) to the list
        results.append(ko_results.drop(columns=["status"]))

    # Create a cocktail condition
    print(f"Running single reaction KOs for cocktail...")
    # All the substrates are available in their ratio in the data
    cocktail_media = minimal_media.copy()
    for index, row in top_10_exometabolites.iterrows():
        cocktail_media[row["exchange_id"]] = (
            TOTAL_UPTAKE
            * row["carbon_concentration"]
            / top_10_exometabolites["carbon_concentration"].sum()
        ) / row["n_c"]
    MODEL.medium = cocktail_media
    # Run the reaction knockouts
    ko_results = cobra.flux_analysis.single_reaction_deletion(MODEL)

    # Set growth to NaN where status is not 'optimal'
    ko_results["growth"] = ko_results.apply(
        lambda row: row["growth"] if row["status"] == "optimal" else np.nan, axis=1
    )

    # Rename the "growth" column to be the name of the substrate
    ko_results = ko_results.rename(columns={"growth": "cocktail"})

    # Convert the 'ids' column from sets to strings
    ko_results["ids"] = ko_results["ids"].apply(lambda x: list(x)[0])

    # Append the results (minus the "status" column) to the list
    results.append(ko_results.drop(columns=["status"]))

    # Combine all dataframes on the 'ids' column.
    # We can do this by setting 'ids' as the index for each df and then concatenating.
    final_df = pd.concat([df.set_index("ids") for df in results], axis=1)

    # Save the final dataframe to a CSV file
    final_df.to_csv(OUT_PATH / "single_reaction_ko_results.csv")


if __name__ == "__main__":
    main()

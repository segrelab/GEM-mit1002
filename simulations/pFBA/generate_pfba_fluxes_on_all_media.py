import json
import os
import pickle

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from gem_utilities import biomass, media

# Define paths relative to the script or project root
# It's better practice to define a project root
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
TESTFILE_DIR = os.path.join(PROJECT_ROOT, "test", "test_files")
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")

# Ensure the results directory exists
os.makedirs(RESULTS_DIR, exist_ok=True)

# Load the media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)

# Load the model
model = cobra.io.read_sbml_model(os.path.join(PROJECT_ROOT, "model.xml"))

# Load the TSV of the growth phenotypes
growth_phenotypes = pd.read_csv(
    os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"),
    sep="\t",
    converters={"met_id": lambda x: x.split(",")},
)

# Loop through the growth phenotpes, and add the carbon source to the
# minimal media, run pFBA and save the fluxes
for index, row in growth_phenotypes.iterrows():
    # Create a safe filename for the media
    media_name = row["minimal_media"] + "_" + row["c_source"]
    file_name = os.path.join(RESULTS_DIR, f"{media_name}_pfba_fluxes.json")
    # If the file already exists, delete it
    if os.path.exists(file_name):
        os.remove(file_name)
    # For debugging, print the media being processed
    print(
        f"Processing media: {row['minimal_media']} with carbon source: {row['c_source']}"
    )
    # If the model cannot grow on the media, skip it
    if row["growth"] == "No":
        continue
    # Get the minimal media
    minimal_media = media_definitions[row["minimal_media"]].copy()
    # Check if the model has an exchange reaction for the metabolite
    if all(
        "EX_" + met_id + "_e0" in [r.id for r in model.reactions]
        for met_id in row["met_id"]
    ):
        # If it does, add the exchange reaction to the minimal media used
        for met_id in row["met_id"]:
            met_obj = model.metabolites.get_by_id(met_id + "_c0")
            # If the metabolite has carbon in it, handle the bounds
            if "C" in met_obj.elements:
                # Set the uptake rate to the equivalent of -10 for glucose
                minimal_media["EX_" + met_id + "_e0"] = 60 / met_obj.elements["C"]
            # If not (e.g. it's a nitrogen source), add an unlimited amount
            else:
                minimal_media["EX_" + met_id + "_e0"] = 1000.0
    # Set the media
    model.medium = media.clean_media(model, minimal_media)
    # Try to run pFBA
    try:
        pfba_solution = cobra.flux_analysis.pfba(model)
    except Exception as e:
        print(f"Error running pFBA for media {row['minimal_media']} with carbon source {row['c_source']}: {e}")
        continue
    # Save the fluxes to a json file for use with Escher
    with open(file_name, "w") as f:
        json.dump(pfba_solution.fluxes.to_dict(), f)

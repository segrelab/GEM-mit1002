import json
import os

import pandas as pd

FILE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(FILE_DIR, "..", "..", "..", "data")

# Load the media database (without carbon sources)
media_db = pd.read_csv(os.path.join(FILE_DIR, "no_c_media_database.tsv"), sep="\t", header=0)

# Load the growth phenotype data
growth_phenotype = pd.read_csv(os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"), sep="\t", header=0)

# Split the media database into the two different media (value in the "medium" column)
media_db_mbm = media_db[media_db["medium"] == "mbm"]
media_db_l1 = media_db[media_db["medium"] == "l1"]

# Need to load in the ModelSEED database first
modelseed_db = json.load(
    open(
        "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Biochemistry/compounds.json"
    )
)
# Convert to a dictionary with the ModelSEED IDs as the keys for easier searching
modelseed_db = {met["id"]: met for met in modelseed_db}


# Write a function to conver the alias strings to a dictionary
def convert_aliases_to_dict(alias_string):
    return {
        alias.split(":")[0]: [ak.strip() for ak in alias.split(":")[1].split(";")]
        for alias in alias_string
        if alias
    }


# Make a new media database with the same columns as the original media database
media_db_c = pd.DataFrame(columns=media_db.columns)

# Loop through the growth phenotype data
for index, row in growth_phenotype.iterrows():
    # Make a copy of the relevant base media
    if row["minimal_media"] == "mbm_media":
        media = media_db_mbm.copy()
    elif row["minimal_media"] == "l1_media":
        media = media_db_l1.copy()
    # Get the medium name and description from the media database
    medium_name = media["medium"].unique()[0]
    medium_description = media["description"].unique()[0]
    # Get the BiGG ID of the carbon source
    aliases = convert_aliases_to_dict(modelseed_db[row["met_id"]]["aliases"])
    if "BiGG" not in aliases:
        print(f"No BiGG ID for {row['met_id']}")
        continue
    else:
        bigg_id = aliases["BiGG"]
        if len(bigg_id) == 0:
            print(f"No BiGG ID for {row['met_id']}")
            continue
        if len(bigg_id) > 1:
            print(f"Multiple BiGG IDs for {row['met_id']}: {bigg_id}")
        bigg_to_use = bigg_id[0]
    # Make a new medium name and description with the BiGG ID of the carbon source
    new_name = medium_name + "_" + bigg_to_use
    new_description = medium_description + " with " + row["c_source"].strip()
    # Add a new row to the media for the carbon source
    new_row_data = {
        "medium": new_name,
        "description": new_description,
        "compound": bigg_to_use,
        "name": row["c_source"],
    }
    new_row = pd.DataFrame(new_row_data, index=[0])
    # Add the new row
    media = pd.concat([media, new_row], ignore_index=True)
    # Set the new medium name and description for all the rows
    media["medium"] = new_name
    media["description"] = new_description
    # Add the new to the new media database
    media_db_c = pd.concat([media_db_c, media])

# Save the new media database
media_db_c.to_csv(os.path.join(FILE_DIR, "media_database.tsv"), sep="\t", index=False)

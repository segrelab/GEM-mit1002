import os

import pandas as pd

# TODO: Set a variable for the one kbase workspace where all my media are stored
media_workspace = "hgsco:narrative_1740349496319"

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(FILE_DIR)), "data")
media = {
    "mbm": {"kbase_id": "mbm.media", "kbase_ws": media_workspace},
    "l1": {"kbase_id": "l1.media", "kbase_ws": media_workspace},
    "bashir_c_free": {"kbase_id": "bashir_c_free.media", "kbase_ws": media_workspace},
    "marine_broth_wo_yeast_and_peptone": {
        "kbase_id": "marine_broth_wo_yeast_and_peptone.media",
        "kbase_ws": media_workspace,
    },
    "promm": {"kbase_id": "promm.media", "kbase_ws": media_workspace},
    "hmb": {"kbase_id": "hmb.media", "kbase_ws": media_workspace},
    "mmb": {"kbase_id": "mmb.media", "kbase_ws": media_workspace},
    "pro99": {"kbase_id": "pro99.media", "kbase_ws": media_workspace},
}

# Load the known growth phenotype data as a pandas DataFrame
growth_phenotype_data = pd.read_csv(
    os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"), sep="\t"
)

# Add a "geneko" ccolumn to the DataFrame and set it to "none"
growth_phenotype_data["geneko"] = "none"

# Add a mediaws column to the DataFrame with the value of "kbase_ws" for
# the media
growth_phenotype_data["mediaws"] = growth_phenotype_data["minimal_media"].apply(
    lambda x: media[x]["kbase_ws"]
)

# Add a media column to the DataFrame with the value of "kbase_id" for
# the media
growth_phenotype_data["media"] = growth_phenotype_data["minimal_media"].apply(
    lambda x: media[x]["kbase_id"]
)

# Change the title of the "met_id" column to "addtlCpd"
growth_phenotype_data.rename(columns={"met_id": "addtlCpd"}, inplace=True)

# Change the values of the growth column from "Yes", "No", and "Unsure"
# to 1 and 0
growth_phenotype_data["growth"] = growth_phenotype_data["growth"].replace(
    {"Yes": 1, "No": 0, "Unsure": 1}
)

# Save just the columns we need to a new DataFrame
growth_phenotype_data = growth_phenotype_data[
    ["geneko", "mediaws", "media", "addtlCpd", "growth"]
]

# Save the DataFrame to a tab-separated file
growth_phenotype_data.to_csv(
    os.path.join(DATA_DIR, "kbase_phenotype_set.tsv"), sep="\t", index=False
)

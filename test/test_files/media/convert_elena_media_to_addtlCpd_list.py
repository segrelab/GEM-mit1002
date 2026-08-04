# I have KBase media definitions for some of Elena's media from Osnat
# But I want to get them as a list of additional compounds to add to the media
# for the phenotype sets

import os

import pandas as pd
from tools.media import MEDIA

# Load the media definitions
media = MEDIA

# Load the KBase media TSV from Osnat's narrative
# Change this whenever you want to run on a new medium
# media_tsv = "/Users/helenscott/Downloads/HMBaa_media 2/HMBaa_media.tsv"
# media_tsv = "/Users/helenscott/Downloads/HMBoligo_media/HMBoligo_media.tsv"
media_tsv = "/Users/helenscott/Downloads/HMBamisug_media/HMBamisug_media.tsv"
media_df = pd.read_csv(media_tsv, sep="\t")

# Get the list of all the comoounds in the media
media_compounds = media_df["compounds"].tolist()

# Subset the list of media compounds to only those that are not in the
# corresponding media definition (e.g. "hmb" for "hmb_aa")
addtl_compounds = [
    c for c in media_compounds if "EX_" + c + "_e0" not in media["hmb"]
]

print(len(addtl_compounds))

# Print the list of additional compounds with semi-colons separating them
print(";".join(addtl_compounds))

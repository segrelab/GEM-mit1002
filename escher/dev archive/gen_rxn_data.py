import os

import cobra

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(FILE_DIR)

import sys

sys.path.insert(0, str(REPO_DIR))

from tools.media import MEDIA  # noqa: E402

# Load the media defintions from the pickle file
media_definitions = MEDIA

# Load the model
model = cobra.io.read_sbml_model(os.path.join(REPO_DIR, "model.xml"))

####################################
# Glucose as sole carbon source
####################################
# Set the medium to minimal medium
model.medium = media_definitions["minimal_glucose"]
# Run the simulation
sol = model.optimize()
# Save the results as a dictionary with the reaction ID as the key and the
# flux value as the value
glc_rxn_data = {}

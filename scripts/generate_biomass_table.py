import os

import cobra
from gem_utilities.biomass import save_biomass_composition_work_table

# Define paths relative to the script or project root
# It's better practice to define a project root
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
RESULTS_DIR = os.path.join(PROJECT_ROOT, "scripts", "results")

# Load the model
model = cobra.io.read_sbml_model(os.path.join(PROJECT_ROOT, "model.xml"))

# Save the biomass composition table
save_biomass_composition_work_table(
    model=model, mets_to_ignore=["cpd11416_c0"], out_dir=RESULTS_DIR
)

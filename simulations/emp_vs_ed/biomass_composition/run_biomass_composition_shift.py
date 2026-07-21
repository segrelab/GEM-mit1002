import os
import pickle as pkl
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from gem_utilities import biomass

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Set the output directory
OUT_DIR = os.path.dirname(os.path.realpath(__file__))


# Load the model
model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

# Load the media definitions
with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
    media_defs = pkl.load(f)

# Set the medium to be the minimal medium with glucose
# FIXME: In the minimal medium, oxygen is 20, so that is limiting, is that ok?
model.medium = media_defs["minimal_glucose"]

# Define key reaction IDs
BIOMASS_RXN = "bio1_biomass"
CO2_EX_RXN = "EX_cpd00011_e0"
GLC_EX_RXN = "EX_cpd00027_e0"
# Define a reaction to indicate flux through the ED pathway
ed_rxn_id = "rxn01477_c0"
# Define a reaction to indicate flux through the EMP
emp_rxn_id = "rxn00558_c0"
# Number of carbon atoms in glucose (to get the uptake in the carbon atom flux)
GLC_N_C = 6

# Define the amount of carbon in biomass
# FIXME: Can't use the N_C_BIOMASS here, becuase it changes with the biomass composition
# From scripts/results/iHS4156_biomass_composition_work_table.csv
N_C_BIOMASS = 42.948  # mmol C

# Define perturbation levels in log space (e.g., ±log(2) for doubling/halving, ±1 log(10) for orders of magnitude)
# Use an odd number of levels to include the base value
log_perturbation_levels = np.linspace(-1, 1, 101)  # ±1 in log10 scale (1/10 to 10x)

# Define metabolites to skip from the biomass composition
components_to_skip = [
    "cpd00001_c0",  # H2O
    "cpd00002_c0",  # ATP
]

# Get the biomass components dictionary
bio_components = model.reactions.get_by_id(BIOMASS_RXN).metabolites

# Unlump the biomass reaction, add that to the model, and remove the old one
# So that there are no lumped components, like "protein" in the biomass reaction
unlumped_biomass = biomass.unlump_biomass(bio_components, model)
# Remove the biomass reaction's old metabolite dictionary
model.reactions.get_by_id(BIOMASS_RXN).subtract_metabolites(
    model.reactions.get_by_id(BIOMASS_RXN).metabolites
)
# Use that dictionary to replace the old biomass reaction
model.reactions.get_by_id(BIOMASS_RXN).add_metabolites(unlumped_biomass, combine=False)

# Prepare the list of biomass components to process
components_to_run = [
    (component, s_coeff)
    for component, s_coeff in unlumped_biomass.items()
    if s_coeff < 0 and component.id not in components_to_skip
]

# Create a MultiIndex for the columns
columns = pd.MultiIndex.from_product(
    [
        [met.name for (met, coeff) in components_to_run],
        ["perturbed_s_coeff", "growth_rate", "cue", "bge", "ed_flux", "emp_flux"],
    ],
    names=["component", "metric"],
)

# Create an empty DataFrame with MultiIndex columns and the perturbation levels as the index
results = pd.DataFrame(index=log_perturbation_levels, columns=columns)


for component_index, (component, s_coeff) in enumerate(components_to_run, start=1):
    print(
        f"Running component {component_index}/{len(components_to_run)}: '{component.name}'"
    )
    print("Progress: ", end="", flush=True)

    # Loop through the multiplier to use
    for run_index, log_level in enumerate(log_perturbation_levels, start=1):
        # Calculate the new stoichiometric coefficient to use
        # Convert to log space, need to use abs() to avoid negative values
        log_s_coeff = np.log10(abs(s_coeff))
        # Apply the perturbation
        log_perturbed_s_coeff = log_s_coeff + log_level
        # Convert back to linear scale, and make negative
        perturbed_s_coeff = -1 * (10**log_perturbed_s_coeff)
        # Set the new stoichiometric coefficient
        model.reactions.get_by_id(BIOMASS_RXN).add_metabolites(
            {component: perturbed_s_coeff}, combine=False
        )

        # Rebalance the stoichiometric coefficicents so the weight is still 1
        # Get the weight with the new coefficient
        new_weight = biomass.calculate_biomass_weight(
            model,
            BIOMASS_RXN,
            mets_to_ignore=["cpd11416_c0"],
            lumped_biomass_components=None,
        )
        if new_weight != 1.000:
            biomass_reaction = model.reactions.get_by_id(BIOMASS_RXN)
            biomass_reaction.add_metabolites(
                {
                    met: coef * (1 / new_weight)
                    for met, coef in biomass_reaction.metabolites.items()
                },
                combine=False,
            )

        # Run pFBA
        sol = cobra.flux_analysis.pfba(model)
        # If the status is infeasible, skip extracting the results
        if sol.status != "optimal":
            growth = None
            cue = None
            gge = None
            bge = None
        else:
            # Extract the growth rate
            growth = sol.fluxes[BIOMASS_RXN]
            # Extract the CO2 release
            co2 = sol.fluxes.get(CO2_EX_RXN, 0.0)

            # Extract the fluxes through the ED and EMP pathways
            ed_flux = sol.fluxes.get(ed_rxn_id, 0.0)
            emp_flux = sol.fluxes.get(emp_rxn_id, 0.0)

            # Get the uptake flux for the exchange ID
            # Absolute value since uptake is negative
            uptake = abs(sol.fluxes.get(GLC_EX_RXN, 0.0))
            # Convert the uptake to mmol C / gDW / h
            uptake_c = uptake * GLC_N_C

            # Calculate the CUE/BGE
            # TODO: Use the helper function
            cue = 1.0 - (co2 / uptake_c)
            # Calculate the BGE (Bacterial Growth Efficiency)
            # FIXME: Can't use the N_C_BIOMASS here, becuase it changes with the biomass composition
            # TODO: Use the helper function
            bge = (N_C_BIOMASS * growth) / ((N_C_BIOMASS * growth) + co2)

        # Save the results to the dataframe
        results.loc[log_level, (component.name, "perturbed_s_coeff")] = (
            perturbed_s_coeff
        )
        results.loc[log_level, (str(component.name), "growth_rate")] = growth
        results.loc[log_level, (str(component.name), "cue")] = cue
        results.loc[log_level, (str(component.name), "bge")] = bge
        results.loc[log_level, (str(component.name), "ed_flux")] = ed_flux
        results.loc[log_level, (str(component.name), "emp_flux")] = emp_flux

        if run_index % 10 == 0 or run_index == len(log_perturbation_levels):
            print("X", end="", flush=True)

    print(f" {run_index}/{len(log_perturbation_levels)} done")

    # Reset the stoichiometric coefficient to the original value
    model.reactions.get_by_id(BIOMASS_RXN).add_metabolites(
        {component: s_coeff}, combine=False
    )

# Save the results as a csv
results.to_csv(OUT_PATH / "results.csv")

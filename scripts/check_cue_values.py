import os
import pickle
import sys
import warnings

import cobra
import matplotlib.pyplot as plt
import pandas as pd
from gem2cue import utils
from gem_utilities import media

# Import the plot styles (has global variables for colors)
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from plot_styles import *

# Set the output directory
FILE_DIR = os.path.dirname(os.path.realpath(__file__))
# Set a folder for the results
OUT_DIR = os.path.join(FILE_DIR, "results")
# Check if the folder exists, if not, create it
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

# ============================================================================
# Load model and define media conditions
# ============================================================================

# Load in the ALT model using COBRApy
alt_cobra = cobra.io.read_sbml_model("model.xml")

# Base minimal media (no carbon source, no O2 — added per condition below)
BASE_MEDIUM = {
    "EX_cpd00067_e0": 1000,  # H+  This has a really big effect on the exudation
    "EX_cpd00058_e0": 1000,  # Cu2+
    "EX_cpd00971_e0": 1000,  # Na+
    "EX_cpd00063_e0": 1000,  # Ca2+
    "EX_cpd00048_e0": 1000,  # Sulfate
    "EX_cpd10516_e0": 1000,  # Fe3+
    "EX_cpd00254_e0": 1000,  # Mg
    "EX_cpd00009_e0": 1000,  # Phosphate
    "EX_cpd00205_e0": 1000,  # K+
    "EX_cpd00013_e0": 1000,  # NH3
    "EX_cpd00099_e0": 1000,  # Cl-
    "EX_cpd00030_e0": 1000,  # Mn2+
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00635_e0": 1000,  # Cbl
    "EX_cpd00034_e0": 1000,  # Zn2+
    "EX_cpd00149_e0": 1000,  # Co2+
}

# Carbon source conditions
# TODO: Use uptake rates based on NMR data
# FIXME: Mix conditions need equivalent carbon amounts
CARBON_SOURCES = {
    "glc": {"EX_cpd00027_e0": 10},  # D-Glucose only
    "ace": {"EX_cpd00029_e0": 30},  # Acetate only
    "glc_heavy_mix": {
        "EX_cpd00027_e0": 6.667,  # 2/3 glucose,
        "EX_cpd00029_e0": 10,
    },  # 1/3 acetate
    "ace_heavy_mix": {
        "EX_cpd00027_e0": 3.333,  # 1/3 glucose,
        "EX_cpd00029_e0": 20,
    },  # 2/3 acetate
}

# Oxygen conditions
O2_CONDITIONS = {
    "inf_o2": {"EX_cpd00007_e0": 1000},  # Unconstrained O2
    "real_o2": {"EX_cpd00007_e0": 20},  # Realistic O2 uptake
}

# Build all media combinations by crossing carbon sources with O2 conditions
media_dict = {
    f"{c_name}_medium_{o2_name}": {**BASE_MEDIUM, **c_fluxes, **o2_fluxes}
    for c_name, c_fluxes in CARBON_SOURCES.items()
    for o2_name, o2_fluxes in O2_CONDITIONS.items()
}

# ============================================================================
# Run FBA and pFBA simulations
# ============================================================================
cobra_results = {}
for name, medium in media_dict.items():
    # Run FBA
    alt_cobra.medium = media.clean_media(alt_cobra, medium)
    fba_result = alt_cobra.optimize()
    cobra_results[name + "_fba"] = fba_result
    # Run pFBA
    pfba_result = cobra.flux_analysis.pfba(alt_cobra)
    cobra_results[name + "_pfba"] = pfba_result

# ============================================================================
# Calculate CUE and other metrics
# ============================================================================
# Load the model and get the exchange reactions
c_ex_rxns = utils.get_c_ex_rxns(alt_cobra)

# Extract the carbon fate results from the FBA results and save them in a DataFrame
results_list = []
for key, fba_result in cobra_results.items():
    # Extract the carbon fates for the solution (both normalized and not normalized)
    c_fates = utils.extract_c_fates_from_solution(
        fba_result, c_ex_rxns, co2_ex_rxn="EX_cpd00011_e0", norm=False
    )
    uptake = c_fates[0]
    co2 = c_fates[1]
    organic_c = c_fates[2]
    biomass = c_fates[3]

    c_fates_norm = utils.extract_c_fates_from_solution(
        fba_result, c_ex_rxns, co2_ex_rxn="EX_cpd00011_e0", norm=True
    )
    co2_norm = c_fates_norm[0]
    organic_c_norm = c_fates_norm[1]
    biomass_norm = c_fates_norm[2]

    # Calculate CUE from the c fates (not using my function)
    cue = 1 - co2 / uptake

    # Calculate GGE from the c fates (not using my function)
    gge = 1 - (co2 + organic_c) / uptake

    # Calculat BGE from the c fates (not using my function)
    bge = biomass / (biomass + co2)

    results_list.append(
        {
            "sim_name": key,
            "oxygen_flux": fba_result.fluxes["EX_cpd00007_e0"],
            "uptake": round(abs(uptake), 3),
            "uptake_norm": 1,  # This is always 1 because the uptake is the reference
            "co2": round(abs(co2), 3),
            "co2_norm": round(abs(co2_norm), 3),
            "organic_c": round(abs(organic_c), 3),
            "organic_c_norm": round(abs(organic_c_norm), 3),
            "biomass": round(abs(biomass), 3),
            "biomass_norm": round(abs(biomass_norm), 3),
            "cue": round(cue, 3),
            "gge": round(gge, 3),
            "bge": round(bge, 3),
        }
    )

# Convert the results to a DataFrame
results = pd.DataFrame(results_list)

# Save the results
results.to_csv(os.path.join(OUT_DIR, "glc_ace_mixes_cue.csv"), index=False)

# ============================================================================
# Generate plots
# ============================================================================
# Stacked bar plot of the carbon fates for the different conditions
data = results.set_index("sim_name")[["biomass", "organic_c", "co2"]]
g = carbon_fates_bar(data)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "glc_ace_mixes_carbon_fates.png"))

# Stacked bar plot of the normalized carbon fates for the different conditions
data_norm = results.set_index("sim_name")[
    ["biomass_norm", "organic_c_norm", "co2_norm"]
]
# Rename the columns to match the function
data_norm.columns = ["biomass", "organic_c", "co2"]
g_norm = carbon_fates_bar(data_norm)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "glc_ace_mixes_carbon_fates_norm.png"))

# Subset the results to only include the pFBA on realistic O2
clean_data = data_norm[data_norm.index.str.contains("real_o2_pfba")]
g_clean = carbon_fates_bar(clean_data)
# Relabel the x tick labels
g_clean.set_xticklabels(
    ["Glucose", "Acetate", "Heavy Glucose Mix", "Heavy Acetate Mix"]
)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "glc_ace_mixes_carbon_fates_norm_real_o2_only.png"))

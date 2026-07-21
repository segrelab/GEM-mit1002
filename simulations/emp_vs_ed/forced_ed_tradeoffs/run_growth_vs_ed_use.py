"""Force increasing flux through the ED pathway on glucose and record the cost.

For a range of forced ED-pathway fluxes (0 -> full glycolytic flux), at both a
limiting O2 supply (20) and excess O2 (1000), run pFBA on minimal glucose medium
and record growth rate, carbon use efficiency (CUE), bacterial growth efficiency
(BGE), and acetate secretion. The full flux vector for each point is also saved
to results/fluxes/ (used for the escher maps).

Converted from plot_growth_vs_ed_use.ipynb; plotting lives in
plot_growth_vs_ed_use.py.
"""

import json
import pickle as pkl
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)
FLUX_PATH = OUT_PATH / "fluxes"
FLUX_PATH.mkdir(exist_ok=True)

# Key reaction IDs (consistent with the other emp_vs_ed scripts)
ED_RXN_ID = "rxn01477_c0"     # KDPG aldolase (eda) -- ED marker
EMP_RXN_ID = "rxn00558_c0"    # PFK -- EMP marker
BIOMASS_RXN_ID = "bio1_biomass"
O2_EX_RXN = "EX_cpd00007_e0"
CO2_EX_RXN = "EX_cpd00011_e0"
GLC_EX_RXN = "EX_cpd00027_e0"
ACETATE_EX_RXN = "EX_cpd00029_e0"

# Carbon atoms per glucose (for CUE)
GLC_N_C = 6
# Carbon atoms in biomass (for BGE)
# From biomass/iHS4156_biomass_composition_work_table.csv, generated with:
# biomass.save_biomass_composition_work_table(model, mets_to_ignore=["cpd11416_c0"])
BIOMASS_C = 42.95

# Forced ED-pathway flux levels and O2 uptake levels to sweep
ED_USE_LEVELS = np.linspace(0, 10, 11)
O2_LEVELS = [20, 1000]


def main():
    # Load the model and set the minimal glucose medium
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
    with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
        media_definitions = pkl.load(f)
    model.medium = media_definitions["minimal_glucose"]

    ed_rxn = model.reactions.get_by_id(ED_RXN_ID)
    emp_rxn = model.reactions.get_by_id(EMP_RXN_ID)

    results = []
    for o2 in O2_LEVELS:
        # Set the O2 uptake rate (negative lower bound = uptake)
        model.reactions.get_by_id(O2_EX_RXN).lower_bound = -o2
        for use in ED_USE_LEVELS:
            # Force flux through the ED pathway via the ED reaction's lower bound
            ed_rxn.lower_bound = use
            # Run pFBA
            sol = cobra.flux_analysis.pfba(model)

            ed_flux = sol.fluxes[ED_RXN_ID]
            emp_flux = sol.fluxes[EMP_RXN_ID]
            growth = sol.fluxes[BIOMASS_RXN_ID]
            co2 = sol.fluxes[CO2_EX_RXN]
            glc_uptake = abs(sol.fluxes[GLC_EX_RXN])

            results.append(
                {
                    "o2": o2,
                    "ed_use": use,
                    "ed_percent": ed_flux / (ed_flux + emp_flux) * 100,
                    "growth_rate": growth,
                    "cue": 1 - (co2 / (glc_uptake * GLC_N_C)),
                    "bge": (BIOMASS_C * growth) / (BIOMASS_C * growth + co2),
                    "acetate_secretion": sol.fluxes.get(ACETATE_EX_RXN, 0.0),
                }
            )

            # Save the full flux vector for this point
            with open(FLUX_PATH / f"fluxes_ed_{use}_o2_{o2}.json", "w") as f:
                json.dump(sol.fluxes.to_dict(), f)

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_PATH / "growth_vs_ed_use_results.csv", index=False)
    print(f"Saved {len(results_df)} rows to {OUT_PATH / 'growth_vs_ed_use_results.csv'}")


if __name__ == "__main__":
    main()

"""The oxygen cost of growing on phenol.

Phenol degradation uses O2 twice as a *co-substrate* (phenol hydroxylase and
catechol 2,3-dioxygenase ring cleavage) before any energy is recovered, on top
of the O2 used as a terminal electron acceptor for respiration. Sugars and
acetate use O2 only for respiration. This script quantifies both effects:

  1. Growth rate vs O2 supply for phenol, glucose, and acetate at a fixed carbon
     uptake -- phenol should be obligately aerobic and steeply O2-limited.
  2. The O2 budget at non-limiting O2, split into "oxygenase" (catabolic
     co-substrate) and "respiratory" (terminal oxidase) O2 per substrate.

Plotting lives in plot_phenol_o2.py.
"""

from pathlib import Path

import cobra
import numpy as np
import pandas as pd

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]

import sys

sys.path.insert(0, str(REPO_ROOT))

from tools.media import MEDIA  # noqa: E402
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

BIOMASS_RXN = "bio1_biomass"
O2_EX_RXN = "EX_cpd00007_e0"

# Substrates to compare (name -> entry in the growth_and_cue substrate panel)
SUBSTRATES = ["Phenol", "Glucose", "Acetate"]

# Fix total carbon uptake so substrates are compared per carbon atom
TOTAL_C_UPTAKE = 60  # mmol C / gDW / hr

# O2 supply levels to sweep (mmol O2 / gDW / hr)
O2_LEVELS = np.linspace(0, 150, 31)
# A non-limiting O2 supply for the budget decomposition
O2_NONLIMITING = 1000

# O2 used as a catabolic co-substrate (oxygenases); coeff = O2 consumed per flux
OXYGENASE_RXNS = {
    "rxn39978_c0": 1.0,  # phenol hydroxylase
    "rxn39979_c0": 1.0,  # phenol hydroxylase
    "rxn00587_c0": 1.0,  # catechol 2,3-dioxygenase (ring cleavage)
}
# O2 used as terminal electron acceptor (respiratory oxidases)
RESPIRATORY_RXNS = {
    "rxn13688_c0": 0.5,  # cytochrome-c oxidase
    "rxn14419_c0": 0.5,  # cytochrome-bd oxidase
    "rxn14422_c0": 0.5,  # cytochrome-bo oxidase
    "rxn14426_c0": 0.5,  # cytochrome oxidase cbb3
}


def build_medium(model, base, exchange_id, uptake, o2):
    med = dict(base)
    med[exchange_id] = uptake
    med[O2_EX_RXN] = o2
    rxn_ids = {r.id for r in model.reactions}
    return {k: v for k, v in med.items() if k in rxn_ids}


def main():
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
    media_defs = MEDIA
    base = media_defs["minimal"]  # minimal medium, no carbon source

    panel = pd.read_csv(
        REPO_ROOT / "simulations" / "growth_and_cue" / "results" / "substrate_panel.csv"
    ).set_index("name")

    # --- 1. Growth vs O2 supply ---
    growth_records = []
    for name in SUBSTRATES:
        row = panel.loc[name]
        uptake = TOTAL_C_UPTAKE / row["n_c"]
        for o2 in O2_LEVELS:
            med = build_medium(model, base, row["exchange_id"], uptake, o2)
            with model:
                model.medium = med
                g = model.slim_optimize()
            growth_records.append(
                {
                    "substrate": name,
                    "o2_supply": o2,
                    "growth_rate": g if (g is not None and not np.isnan(g)) else 0.0,
                }
            )
    pd.DataFrame(growth_records).to_csv(OUT_PATH / "growth_vs_o2.csv", index=False)

    # --- 2. O2 budget at non-limiting O2 ---
    budget_records = []
    for name in SUBSTRATES:
        row = panel.loc[name]
        uptake = TOTAL_C_UPTAKE / row["n_c"]
        med = build_medium(model, base, row["exchange_id"], uptake, O2_NONLIMITING)
        with model:
            model.medium = med
            sol = cobra.flux_analysis.pfba(model)
            f = sol.fluxes
            oxygenase_o2 = sum(f.get(r, 0.0) * c for r, c in OXYGENASE_RXNS.items())
            respiratory_o2 = sum(f.get(r, 0.0) * c for r, c in RESPIRATORY_RXNS.items())
            total_o2 = abs(f.get(O2_EX_RXN, 0.0))
            budget_records.append(
                {
                    "substrate": name,
                    "growth_rate": f[BIOMASS_RXN],
                    "total_o2": total_o2,
                    "oxygenase_o2": oxygenase_o2,
                    "respiratory_o2": respiratory_o2,
                    "other_o2": max(total_o2 - oxygenase_o2 - respiratory_o2, 0.0),
                    "total_o2_per_c": total_o2 / TOTAL_C_UPTAKE,
                    "oxygenase_o2_per_c": oxygenase_o2 / TOTAL_C_UPTAKE,
                    "respiratory_o2_per_c": respiratory_o2 / TOTAL_C_UPTAKE,
                }
            )
    pd.DataFrame(budget_records).to_csv(OUT_PATH / "o2_budget.csv", index=False)
    print("Saved growth_vs_o2.csv and o2_budget.csv to", OUT_PATH)


if __name__ == "__main__":
    main()

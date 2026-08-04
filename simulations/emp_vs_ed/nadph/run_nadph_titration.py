"""Titrate NADPH demand and watch how the cell satisfies it.

We add an O2-coupled "redox maintenance" reaction that consumes NADPH
(NADPH + H+ + 0.5 O2 -> NADP+ + H2O, fully mass/charge balanced) and force an
increasing amount of flux through it. For each forced level we run pFBA
(maximizing biomass) and record how much of the demand is met by the ED
pathway versus the competing NADPH sources (oxidative PPP, NADP-isocitrate
dehydrogenase, malic enzyme, transhydrogenase).

The point: on a normally EMP-preferring substrate (glucose), does rising redox
demand alone pull flux into ED? An ED-obligate uronic acid (galacturonate) is
included as a control, where ED is already required regardless of NADPH demand.
"""

from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from gem_utilities import media as media_utils

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]

import sys

sys.path.insert(0, str(REPO_ROOT))

from tools.media import MEDIA  # noqa: E402
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# --- Key reaction IDs (kept consistent with the other emp_vs_ed scripts) ---
# Marker for flux through the ED pathway (KDPG aldolase, eda)
ED_RXN_ID = "rxn01477_c0"
# Marker for flux through the EMP pathway (phosphofructokinase)
EMP_RXN_ID = "rxn00558_c0"
# Biomass reaction
BIOMASS_RXN_ID = "bio1_biomass"
# CO2 exchange (for CUE)
CO2_EX_RXN = "EX_cpd00011_e0"

# Competing NADPH-producing reactions to track across the sweep.
# In this model ED draws its NADPH from the shared G6PDH step (rxn00604) and
# then diverges from the oxidative PPP at 6-phosphogluconate; the decarboxylating
# 6-phosphogluconate dehydrogenase (rxn01115) is the oxPPP branch that loses a
# carbon as CO2. The others are the classic alternative NADPH valves.
NADPH_SOURCE_RXNS = {
    "G6PDH (rxn00604)": "rxn00604_c0",
    "6PGDH / oxPPP (rxn01115)": "rxn01115_c0",
    "NADP-IDH (rxn01387)": "rxn01387_c0",
    "Malic enzyme (rxn00161)": "rxn00161_c0",
    "Transhydrogenase (rxn09295)": "rxn09295_c0",
}

# --- Metabolite IDs for the NADPH drain reaction ---
NADPH = "cpd00005_c0"
NADP = "cpd00006_c0"
PROTON = "cpd00067_c0"
O2 = "cpd00007_c0"
H2O = "cpd00001_c0"
DRAIN_RXN_ID = "NADPH_redox_maint_c0"

# Oxygen exchange reaction
O2_EX_RXN = "EX_cpd00007_e0"

# Total carbon uptake to fix across substrates (mmol C / gDW / hr), matching the
# forced-routing analysis so growth rates are comparable.
TOTAL_UPTAKE = 60

# O2 uptake bounds to sweep over (matches forced_ed_tradeoffs: limited vs excess)
O2_LEVELS = [20, 1000]

# NADPH drain levels to force (mmol / gDW / hr). Infeasible levels are recorded
# as NaN so the plotting script can simply drop them.
DRAIN_LEVELS = np.linspace(0, 20, 41)


def add_nadph_drain(model):
    """Add the O2-coupled NADPH redox-maintenance reaction to the model."""
    rxn = cobra.Reaction(DRAIN_RXN_ID)
    rxn.name = "NADPH redox maintenance (O2-coupled drain)"
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    rxn.add_metabolites(
        {
            model.metabolites.get_by_id(NADPH): -1.0,
            model.metabolites.get_by_id(PROTON): -1.0,
            model.metabolites.get_by_id(O2): -0.5,
            model.metabolites.get_by_id(NADP): 1.0,
            model.metabolites.get_by_id(H2O): 1.0,
        }
    )
    model.add_reactions([rxn])
    # Sanity check: the reaction should be mass and charge balanced.
    imbalance = rxn.check_mass_balance()
    if imbalance:
        raise ValueError(f"NADPH drain reaction is not balanced: {imbalance}")
    return rxn


def main():
    # Load the model and add the NADPH drain reaction
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
    add_nadph_drain(model)

    # Load the media definitions
    media_defs = MEDIA

    # Load the substrate panel and pick glucose (EMP-preferring) plus
    # galacturonic acid (ED-obligate control)
    substrate_df = pd.read_csv(
        REPO_ROOT / "simulations" / "growth_and_cue" / "results" / "substrate_panel.csv"
    )
    substrates = substrate_df[
        substrate_df["name"].isin(["Glucose", "Galacturonic Acid"])
    ]

    records = []
    for _, row in substrates.iterrows():
        substrate_name = row["name"]
        substrate_uptake = TOTAL_UPTAKE / row["n_c"]

        for o2_level in O2_LEVELS:
            # Build the medium: minimal base + the chosen carbon source, with O2
            # set to the current level.
            media = media_utils.clean_media(model, media_defs["minimal"])
            media[row["exchange_id"]] = substrate_uptake
            media[O2_EX_RXN] = o2_level
            model.medium = media

            for drain_level in DRAIN_LEVELS:
                with model:
                    # Force at least `drain_level` of NADPH oxidation
                    model.reactions.get_by_id(DRAIN_RXN_ID).lower_bound = drain_level

                    try:
                        sol = cobra.flux_analysis.pfba(model)
                        status = sol.status
                    except Exception:
                        status = "infeasible"

                    rec = {
                        "substrate": substrate_name,
                        "o2_level": o2_level,
                        "nadph_drain": drain_level,
                        "status": status,
                    }

                    if status == "optimal":
                        f = sol.fluxes
                        growth = f[BIOMASS_RXN_ID]
                        ed = f.get(ED_RXN_ID, 0.0)
                        emp = f.get(EMP_RXN_ID, 0.0)
                        co2 = f.get(CO2_EX_RXN, 0.0)
                        uptake_c = abs(f.get(row["exchange_id"], 0.0)) * row["n_c"]

                        rec["growth_rate"] = growth
                        rec["ed_flux"] = ed
                        rec["emp_flux"] = emp
                        rec["percent_ed_flux"] = (
                            ed / (ed + emp) if (ed + emp) != 0 else np.nan
                        )
                        rec["co2_flux"] = co2
                        rec["cue"] = (
                            1.0 - (co2 / uptake_c) if uptake_c != 0 else np.nan
                        )
                        rec["drain_flux"] = f.get(DRAIN_RXN_ID, 0.0)
                        for label, rid in NADPH_SOURCE_RXNS.items():
                            rec[label] = f.get(rid, 0.0)
                    else:
                        # Mark unmet columns as NaN
                        for key in [
                            "growth_rate", "ed_flux", "emp_flux", "percent_ed_flux",
                            "co2_flux", "cue", "drain_flux",
                        ] + list(NADPH_SOURCE_RXNS.keys()):
                            rec[key] = np.nan

                    records.append(rec)

        print(f"Done: {substrate_name}")

    results = pd.DataFrame.from_records(records)
    results.to_csv(OUT_PATH / "nadph_titration_results.csv", index=False)
    print(f"Saved {len(results)} rows to {OUT_PATH / 'nadph_titration_results.csv'}")


if __name__ == "__main__":
    main()

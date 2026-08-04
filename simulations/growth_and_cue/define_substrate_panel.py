import json
from pathlib import Path
from typing import Optional

import cobra
import pandas as pd
from gem_utilities import media as media_utils

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]

import sys

sys.path.insert(0, str(REPO_ROOT))

from tools.media import MEDIA  # noqa: E402
DATA_DIR = REPO_ROOT / "data"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Entry point into central metabolism = the first central-metabolism
# intermediate the substrate's catabolism produces.
ENTRY_POINT = {
    "cpd00027": "Glycolysis",  # Glucose      → G6P
    "cpd00108": "Glycolysis",  # Galactose    → G1P → G6P
    "cpd00080": "Glycolysis",  # Glycerol-3-P → DHAP
    "cpd00035": "Pyruvate",  # Alanine      → Pyr
    "cpd00033": "Pyruvate",  # Glycine      → Ser → Pyr
    "cpd23538": "Pyruvate",  # DHPS         → Pyr + sulfite
    "cpd00029": "Acetyl-CoA",  # Acetate      → AcCoA
    "cpd00797": "Acetyl-CoA",  # 3-HB         → 2× AcCoA
    "cpd00107": "Acetyl-CoA",  # Leucine      → AcCoA
    "cpd00039": "Acetyl-CoA",  # Lysine       → AcCoA
    "cpd00023": "TCA — α-KG",  # Glutamate    → α-KG
    "cpd00129": "TCA — α-KG",  # Proline      → Glu → α-KG
    "cpd00051": "TCA — α-KG",  # Arginine     → Glu → α-KG
    "cpd00041": "TCA — OAA",  # Aspartate    → OAA
    "cpd00156": "TCA — Succinyl-CoA",  # Valine       → succinyl-CoA
    "cpd00322": "TCA — Succinyl-CoA",  # Isoleucine   → succinyl-CoA
    "cpd00123": "TCA — Succinyl-CoA",  # KIC          → succinyl-CoA
    "cpd00069": "Aromatic catabolism",  # Tyrosine     → fumarate + AcAcCoA
    "cpd00127": "Aromatic catabolism",  # Phenol       → succinyl-CoA + AcCoA
    "cpd01384": "ED - KDPG",  # D-Mannuronate → KDPG
    "cpd00280": "ED - KDPG",  # D-Galacturonate → KDPG
}

# Define the biomass precurosors and their weights in the pool reaction
PRECURSORS = {
    "cpd00079_c0": 1.0,  # G6P
    "cpd00072_c0": 1.0,  # F6P
    "cpd00101_c0": 1.0,  # R5P
    "cpd00236_c0": 1.0,  # E4P
    "cpd00102_c0": 1.0,  # GAP / triose-P
    "cpd00169_c0": 1.0,  # 3PG
    "cpd00061_c0": 1.0,  # PEP
    "cpd00020_c0": 1.0,  # pyruvate
    "cpd00022_c0": 1.0,  # acetyl-CoA
    "cpd00032_c0": 1.0,  # OAA
    "cpd00024_c0": 1.0,  # alpha-KG
    "cpd00078_c0": 1.0,  # succinyl-CoA
}
COA = "cpd00010_c0"
THIOESTERS = {
    "cpd00022_c0",  # Acetyl-CoA
    "cpd00078_c0",  # succinyl-CoA
}

# Paramaters for that atp_cost function
NTP = [
    "cpd00002_c0",  # ATP
    "cpd00038_c0",  # GTP
]
ATP_SYNTHASE = "rxn08173_c0"
ATPM = "rxn00062_c0"
BIOMASS_RXN = "bio1_biomass"
POOL_RATE = 0.01  # fixed precursor-pool draw
TOTAL_UPTAKE = 60  # mmol C / gDW / hr


def main():
    # Load the model
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    # Load the media definitions
    media_defs = MEDIA

    # Load the known growth phenotypes table with predicted results
    df = pd.read_csv(DATA_DIR / "known_growth_phenotypes.tsv", sep="\t")

    # Filter the growth phenotypes to keep:
    # substrates where growth is observed
    df = df[df["growth"] == "Yes"]
    # substrates that don't have a comma in the met_id (those are cases where multiple substrates were added together)
    df = df[~df["met_id"].astype(str).str.contains(",", na=True)].copy()

    # Drop rows with duplicate met_id values
    df = df.drop_duplicates(subset=["met_id"], keep="first")

    # Drop the unnecessary columns
    # Keep only "c_source" and "met_id"
    df = df.drop(df.columns.difference(["c_source", "met_id"]), axis=1)

    # Add pyruvate to the substrate panel
    # Franzi tested it in conjunction with other substrates, but not alone
    # Add aspartate and glycine to the substrate panel
    # They were measured in the Prochlorococcus exometabolome, but not tested for Amac growth
    # Only add if they aren't already there
    for name, met_id in [
        ("Pyruvate", "cpd00020"),
        ("Aspartate", "cpd00041"),
        ("Glycine", "cpd00033"),
    ]:
        if met_id not in df["met_id"].values:
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        [
                            {
                                "minimal_media": "l1",
                                "c_source": name,
                                "met_id": met_id,
                                "growth": "Yes",
                                "pro_exomet": "Yes",
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

    # Get a list of reaction IDs in the model for later use
    rxn_ids = {r.id for r in model.reactions}

    # Loop through the substrates (rows) and build a record for each one with the information we want to keep
    records = []
    for _, row in df.iterrows():
        # Extract information from the row, stripping whitespace and converting to strings
        met_id = str(row["met_id"]).strip()
        c_source = str(row["c_source"]).strip()
        media_key = str(row["minimal_media"]).strip()

        # Check that the exchange reaction for this substrate exists in the model
        ex_id = f"EX_{met_id}_e0"
        if ex_id not in rxn_ids:
            print(f"  SKIP (no exchange rxn)  : {c_source} ({met_id})")
            continue

        # Get the metabolite
        met = model.metabolites.get_by_id(f"{met_id}_e0")

        # Get basic information about the metabolite
        n_c = met.elements.get("C", 0)
        model_name = met.name[:-5]  # drop " [e0]" suffix

        # Build the media (minimal + substrate)
        media = media_utils.clean_media(model, media_defs["minimal"])
        # Set the substrate uptake to be the total uptake divided by the number of carbons
        # So that every substrate has the same amount of carbon available
        media[ex_id] = TOTAL_UPTAKE / n_c
        # Make the oxygen level unlimited
        media["EX_cpd00007_e0"] = 1000
        # Set the media for the model
        model.medium = media

        # Calculate the NOSC
        nosc = calculate_nosc(met.elements, met.charge)

        # Run pFBA
        sol = cobra.flux_analysis.pfba(model)
        # Get the oxygen uptake rate
        o2_uptake = abs(sol.fluxes.get("EX_cpd00007_e0", 0.0))

        # Calculate the ATP cost of using this substrate
        print(f"Running the ATP Cost calculation for {c_source} ({met_id})...")
        atp_cost = calculate_atp_cost(model, ex_id, pool_rate=POOL_RATE)

        # Add the record for this substrate to the list of records
        records.append(
            {
                "name": c_source,
                "name_in_model": model_name,
                "met_id": met_id,
                "exchange_id": ex_id,
                "n_c": n_c,
                "entry_point": ENTRY_POINT.get(met_id, "Other"),
                "nosc": nosc,
                "o2_saturation": o2_uptake,
                "atp_cost": atp_cost,
            }
        )

    # Convert the list of records to a DataFrame
    substrate_df = pd.DataFrame(records)

    # Save the substrate df
    substrate_df.to_csv(OUT_PATH / "substrate_panel.csv", index=False)


def calculate_nosc(elements: dict, charge: float) -> Optional[float]:
    # Extract the number of atoms of each element, defaulting to 0 if not present
    C = elements.get("C", 0)
    H = elements.get("H", 0)
    N = elements.get("N", 0)
    O = elements.get("O", 0)
    P = elements.get("P", 0)
    S = elements.get("S", 0)

    # Calulate NOSC
    nosc = 4 - ((4 * C) + H - (3 * N) - (2 * O) + (5 * P) - (2 * S) - charge) / C
    return nosc


def net_ntp(model, sol, met_ids, exclude=()):
    """Net high-energy phosphate produced (negative = net consumed),
    excluding reactions in `exclude` (e.g. the ATP synthase)."""
    total = 0.0
    for mid in met_ids:
        m = model.metabolites.get_by_id(mid)
        for r in m.reactions:
            if r.id not in exclude:
                total += sol.fluxes[r.id] * r.metabolites[m]
    return total


def calculate_atp_cost(model, c_source_ex_rxn, pool_rate):
    with model:
        # Minimise substrate uptake to force the most carbon-efficient route
        # Max of a negative uptake flux => least uptake
        model.objective = c_source_ex_rxn
        model.objective_direction = "max"

        # run pFBA
        sol = cobra.flux_analysis.pfba(model)

        # For debugging
        # Save the fluxes as a JSON file
        # fluxes = sol.fluxes.to_dict()
        # with open(OUT_PATH / f"{c_source_ex_rxn}_atp_cost_flux.json", "w") as f:
        #     json.dump(fluxes, f, indent=2)

        # Calculate the net NTP production
        net = net_ntp(model, sol, NTP, exclude={ATP_SYNTHASE})

    return -net / pool_rate  # negative net => consumed => positive cost


if __name__ == "__main__":
    main()

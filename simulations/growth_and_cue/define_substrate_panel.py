import re
import warnings
from pathlib import Path
from typing import Optional
import pickle as pkl

import cobra
from gem_utilities import media as media_utils
import pandas as pd

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

TOTAL_UPTAKE = 60  # mmol C / gDW / hr
BIOMASS_RXN = "bio1_biomass"
CO2_EX_RXN = "EX_cpd00011_e0"

# Exchange metabolites whose max |flux| across substrates is below this are
# lumped into a single grey "Other" segment (keeps the trace-ion colours out).
EX_FLUX_THRESHOLD = 1.0  # mmol / gDW / hr

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
}

# Paramaters for that atp_cost function
NTP = ["cpd00002_c0", "cpd00038_c0"]  # ATP + GTP
ATP_SYNTHASE = "rxn08173_c0"  # VERIFY against your model (see below)
UPTAKE = 10.0  # forced substrate uptake (mmol/gDW/h)
POOL_RATE = 1.0  # fixed precursor-pool draw


def main():
    print("Loading model...")
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    print("\nBuilding substrate panel...")
    substrate_df = load_substrates(model)

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


def atp_cost(model, info, pool_rate=0.1):
    with model:
        configure_medium(
            model, info["ex"], UPTAKE, force_exact=False
        )  # uptake free, we minimise it
        model.reactions.get_by_id(BIOMASS_RXN).bounds = (0, 0)
        pool = cobra.Reaction("DM_pool")
        pool.add_metabolites(
            {model.metabolites.get_by_id(m): -w for m, w in PRECURSORS.items()}
        )
        pool.bounds = (pool_rate, pool_rate)  # force the conversion
        model.add_reactions([pool])
        # carbon-efficient route: minimise substrate uptake
        model.objective = info["ex"]
        model.objective_direction = (
            "max"  # max of a negative uptake flux => least uptake
        )
        sol = cobra.flux_analysis.pfba(model)
        net = net_ntp(model, sol, NTP, exclude={ATP_SYNTHASE})
    return -net / pool_rate  # negative net => consumed => positive cost


def load_substrates(model: cobra.Model) -> pd.DataFrame:
    """Build the single-substrate panel from the known-growth-phenotypes TSV."""
    df = pd.read_csv(TEST_FILE_DIR / "known_growth_phenotypes.tsv", sep="\t")
    df = df[df["growth"] == "Yes"].copy()
    df = df[~df["met_id"].astype(str).str.contains(",", na=True)].copy()

    if "pro_exomet" in df.columns:
        df["pro_exomet"] = df["pro_exomet"].fillna("No")
    else:
        df["pro_exomet"] = "No"

    # Keep the Pro-exometabolite-annotated row when a substrate appears in
    # multiple media contexts (stable sort puts "Yes" first, then dedup).
    df = (
        df.sort_values(
            "pro_exomet",
            key=lambda s: s.map({"Yes": 0}).fillna(1),
            ascending=True,
            kind="stable",
        )
        .drop_duplicates(subset="met_id", keep="first")
        .copy()
    )

    # Aspartate and glycine were measured in the Pro exometabolome; add if absent.
    for name, met_id in [("Aspartate", "cpd00041"), ("Glycine", "cpd00033")]:
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

    rxn_ids = {r.id for r in model.reactions}
    records = []
    for _, row in df.iterrows():
        met_id = str(row["met_id"]).strip()
        c_source = str(row["c_source"]).strip()
        ex_id = f"EX_{met_id}_e0"
        media_key = str(row["minimal_media"]).strip()

        if ex_id not in rxn_ids:
            print(f"  SKIP (no exchange rxn)  : {c_source} ({met_id})")
            continue

        # Get the metabolite
        met = model.metabolites.get_by_id(f"{met_id}_e0")

        # Get information about the metabolite
        n_c = met.elements.get("C", 0)
        model_name = met.name[:-5]  # drop " [e0]" suffix
        nosc = calculate_nosc(met.elements, met.charge)
        atp_cost = atp_cost(model)

        records.append(
            {
                "name": c_source,
                "name_in_model": model_name,
                "met_id": met_id,
                "exchange_id": ex_id,
                "media_key": media_key,
                "n_c": n_c,
                "entry_point": ENTRY_POINT.get(met_id, "Other"),
                "nosc": nosc,
                "atp_cost": atp_cost,
            }
        )

    substrate_df = pd.DataFrame(records)
    return substrate_df


if __name__ == "__main__":
    main()

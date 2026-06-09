"""MIT1002 growth rate and carbon-use efficiency (CUE) across a broad substrate panel.

Runs pFBA on a curated single-substrate panel at fixed total carbon uptake
(60 mmol C / gDW / hr) and reports, per substrate:
  - growth rate (biomass flux)
  - CUE = 1 - CO2_secreted / C_taken_up

Bars are coloured by the substrate's entry point into central metabolism and
sorted by growth rate. Self-contained: re-runs the simulations, no dependence
on the (retired) PCA analysis.
"""

import re
import warnings
from pathlib import Path
from typing import Optional
import pickle as pkl

import cobra
import cobra.flux_analysis
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


def main():
    print("Loading model...")
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    print("Loading media definitions...")
    with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
        media_defs = pkl.load(f)

    print("\nBuilding substrate panel...")
    substrate_df = load_substrates(model, media_defs)

    # Save the substrate df
    substrate_df.to_csv(OUT_PATH / "substrate_panel.csv", index=False)

    print("\nRunning pFBA simulations...")
    summary_df, ex_records = run_pfba(model, media_defs, substrate_df)
    print(f"\nSuccessful: {len(summary_df)}/{len(substrate_df)} substrates")

    # Save the results
    summary_df.to_csv(OUT_PATH / "growth_and_cue.csv")

    # Extract the exchange fluxes
    order = summary_df.sort_values("growth_rate", ascending=False).index.tolist()
    ex_df = build_exchange_df(model, ex_records, order, EX_FLUX_THRESHOLD)
    ex_df.to_csv(OUT_PATH / "exchange_fluxes.csv")


def count_carbons(formula: str) -> Optional[int]:
    """Number of carbon atoms in a molecular formula string."""
    if not formula:
        return None
    m = re.search(r"C(\d*)", formula)
    if m:
        n = m.group(1)
        return int(n) if n else 1
    return 0


def load_substrates(model: cobra.Model, media_defs: dict) -> pd.DataFrame:
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
        if media_key not in media_defs:
            print(f"  SKIP (unknown media '{media_key}'): {c_source}")
            continue

        try:
            n_c = count_carbons(model.metabolites.get_by_id(f"{met_id}_e0").formula)
            model_name = model.metabolites.get_by_id(f"{met_id}_e0").name[
                :-5
            ]  # drop " [e0]" suffix
        except KeyError:
            n_c = None
            for met in model.metabolites:
                if met.id.startswith(met_id):
                    n_c = count_carbons(met.formula)
                    model_name = met.name
                    break
        if not n_c:
            print(f"  SKIP (can't determine n_c): {c_source}")
            continue

        records.append(
            {
                "name": c_source,
                "name_in_model": model_name,
                "met_id": met_id,
                "exchange_id": ex_id,
                "media_key": media_key,
                "n_c": n_c,
                "entry_point": ENTRY_POINT.get(met_id, "Other"),
            }
        )

    substrate_df = pd.DataFrame(records)
    print(f"\nSubstrate panel: {len(substrate_df)} substrates")
    return substrate_df


def run_pfba(model, media_defs, substrate_df) -> tuple:
    """Run pFBA per substrate.

    Returns (summary_df, ex_records) where summary_df holds growth + CUE and
    ex_records maps {substrate: {exchange_rxn_id: flux}} (non-zero fluxes only).
    """
    ex_rxn_ids = [r.id for r in model.reactions if r.id.startswith("EX_")]
    rows = []
    ex_records = {}
    for _, row in substrate_df.iterrows():
        name = row["name"]
        media = media_utils.clean_media(model, media_defs[row["media_key"]])
        media[row["exchange_id"]] = TOTAL_UPTAKE / row["n_c"]
        with model:
            model.medium = media
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sol = cobra.flux_analysis.pfba(model)
                growth = sol.fluxes[BIOMASS_RXN]
                if growth < 1e-6:
                    print(f"  WARN (near-zero growth)  : {name}  mu={growth:.5f}")
                    continue
                co2 = sol.fluxes.get(CO2_EX_RXN, 0.0)
                cue = 1.0 - (co2 / TOTAL_UPTAKE)
                rows.append(
                    {
                        "substrate": name,
                        "met_id": row["met_id"],
                        "growth_rate": growth,
                        "co2_flux": co2,
                        "cue": cue,
                        "entry_point": row["entry_point"],
                    }
                )
                ex_records[name] = {
                    rid: sol.fluxes[rid]
                    for rid in ex_rxn_ids
                    if abs(sol.fluxes[rid]) > 1e-9
                }
                print(f"  OK   {name:28s}  mu={growth:.4f}  CUE={cue:.3f}")
            except Exception as exc:
                print(f"  FAIL {name}: {exc}")
    return pd.DataFrame(rows).set_index("substrate"), ex_records


def build_exchange_df(model, ex_records, substrate_order, threshold):
    """{substrate: {ex_rxn: flux}} -> DataFrame (substrate x metabolite name).

    Columns are renamed from exchange-reaction id to the metabolite name.
    Trace metabolites (max |flux| < threshold) are collapsed into 'Other'.
    """
    df = pd.DataFrame(ex_records).T.reindex(substrate_order).fillna(0.0)
    # Rename the column with the metabolite name instead of the reaction ID
    rename = {}
    for rid in df.columns:
        met = next(iter(model.reactions.get_by_id(rid).metabolites))
        met_name = met.name
        # If the name ends with " [e0]", drop that suffix for cleaner labels (extracellular met)
        suffix = " [e0]"
        if met_name.endswith(suffix):
            met_name = met_name[: -len(suffix)]
        rename[rid] = met_name
    df = df.rename(columns=rename)

    # collapse duplicate metabolite-name columns if any arise (version-safe)
    if df.columns.duplicated().any():
        df = df.T.groupby(level=0).sum().T

    if threshold > 0:
        keep = [c for c in df.columns if df[c].abs().max() >= threshold]
        small = [c for c in df.columns if c not in keep]
        out = df[keep].copy()
        if small:
            out["Other"] = df[small].sum(axis=1)
        df = out
    return df


if __name__ == "__main__":
    main()

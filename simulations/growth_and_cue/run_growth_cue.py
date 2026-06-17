"""MIT1002 growth rate and carbon-use efficiency (CUE) across a broad substrate panel.

Runs pFBA on a curated single-substrate panel at fixed total carbon uptake
(60 mmol C / gDW / hr) and reports, per substrate:
  - growth rate (biomass flux)
  - CUE = 1 - CO2_secreted / C_taken_up

Bars are coloured by the substrate's entry point into central metabolism and
sorted by growth rate. Self-contained: re-runs the simulations, no dependence
on the (retired) PCA analysis.
"""

import pickle as pkl
import re
import warnings
from pathlib import Path
from typing import Optional

import cobra
import cobra.flux_analysis
import pandas as pd
from gem_utilities import media as media_utils

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)
FLUX_PATH = OUT_PATH / "fluxes"
FLUX_PATH.mkdir(exist_ok=True)

# Set the total carbon uptake to use
TOTAL_UPTAKE = 60  # mmol C / gDW / hr

# Exchange metabolites whose max |flux| across substrates is below this are
# lumped into a single grey "Other" segment (keeps the trace-ion colours out).
EX_FLUX_THRESHOLD = 1.0  # mmol / gDW / hr

# Define key reaction IDs
BIOMASS_RXN = "bio1_biomass"
CO2_EX_RXN = "EX_cpd00011_e0"


def main():
    print("Loading model...")
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    print("Loading media definitions...")
    with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
        media_defs = pkl.load(f)

    print("Loading the substrate panel...")
    substrate_df = pd.read_csv(OUT_PATH / "substrate_panel.csv")

    print("\nRunning pFBA simulations...")
    summary_df, ex_records = run_pfba(
        model, media_defs, substrate_df, save_fluxes=True, fluxes_path=FLUX_PATH
    )
    print(f"\nSuccessful: {len(summary_df)}/{len(substrate_df)} substrates")

    # Save the results
    summary_df.to_csv(OUT_PATH / "growth_and_cue.csv")

    # Extract the exchange fluxes
    order = summary_df.sort_values("growth_rate", ascending=False).index.tolist()
    ex_df = build_exchange_df(model, ex_records, order, EX_FLUX_THRESHOLD)
    ex_df.to_csv(OUT_PATH / "exchange_fluxes.csv")


def run_pfba(
    model, media_defs, substrate_df, save_fluxes=False, fluxes_path=None
) -> tuple:
    """Run pFBA per substrate.

    Returns (summary_df, ex_records) where summary_df holds growth + CUE and
    ex_records maps {substrate: {exchange_rxn_id: flux}} (non-zero fluxes only).
    """
    ex_rxn_ids = [r.id for r in model.reactions if r.id.startswith("EX_")]
    rows = []
    ex_records = {}
    for _, row in substrate_df.iterrows():
        name = row["name"]
        media = media_utils.clean_media(model, media_defs["minimal"])
        # Add the carbon source to the media
        media[row["exchange_id"]] = TOTAL_UPTAKE / row["n_c"]
        # Change the oxygen level to be unlimited
        media["EX_cpd00007_e0"] = 1000
        with model:
            model.medium = media
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sol = cobra.flux_analysis.pfba(model)

                # Save all of the fluxes (if requested)
                if save_fluxes:
                    sol.fluxes.to_json(fluxes_path / f"{name}.json")

                # Extract specific fluxes
                growth = sol.fluxes[BIOMASS_RXN]
                if growth < 1e-6:
                    print(f"  WARN (near-zero growth)  : {name}  mu={growth:.5f}")
                    continue
                co2 = sol.fluxes.get(CO2_EX_RXN, 0.0)
                # Get the uptake flux for the exchange ID
                # Absolute value since uptake is negative
                uptake = abs(sol.fluxes.get(row["exchange_id"], 0.0))
                # Convert the uptake to mmol C / gDW / h
                uptake_c = uptake * row["n_c"]

                # Calculate the CUE/BGE
                # TODO: Use the helper function
                cue = 1.0 - (co2 / uptake_c)
                # TODO: Add BGE calculation here

                # Extract the exchange fluxes
                ex_records[name] = {
                    rid: sol.fluxes[rid]
                    for rid in ex_rxn_ids
                    if abs(sol.fluxes[rid]) > 1e-9
                }

                # Add the substrate results to the full results
                rows.append(
                    {
                        "substrate": name,
                        "met_id": row["met_id"],
                        "growth_rate": growth,
                        "co2_flux": co2,
                        "cue": cue,
                    }
                )
                # Print a status message
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

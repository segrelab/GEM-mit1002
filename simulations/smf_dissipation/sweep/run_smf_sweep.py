"""
SMF dissipation sweep: effect of forced Na+ import (flagellar motor load)
on growth and overflow metabolism across carbon sources.

Adds a flagellar Na+ import reaction (cpd00971_e0 -> cpd00971_c0) and sweeps
its lower bound to simulate SMF dissipation. Tracks growth, NaNQR flux, Na+
exchange with medium, and organic acid secretion on a panel of substrates
that differ in whether they use Na+-coupled symporters for uptake.

The key sanity check: does forcing Na+ import actually cost the cell energy
(via increased NaNQR) rather than just flowing trivially from medium? This is
verified by tracking EX_cpd00971_e0 -- if it stays near zero the Na+ is being
cycled internally through NaNQR, which is the real SMF cost.
"""

import json
import warnings
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]

import sys

sys.path.insert(0, str(REPO_ROOT))

from tools.media import MEDIA  # noqa: E402
OUT_PATH = FILE_PATH / "results"
FLUX_DIR = OUT_PATH / "flux_distributions"
OUT_PATH.mkdir(exist_ok=True)
FLUX_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TOTAL_UPTAKE = 60   # mmol C / gDW / hr (matches other analyses in this repo)
BIOMASS_RXN  = "bio1_biomass"
NANQR_RXN    = "ec7211_c0"          # Na+-translocating NADH:quinone oxidoreductase
NA_EX_RXN    = "EX_cpd00971_e0"     # Na+ medium exchange (positive = secretion)
ATPSYN_RXN   = "rxn08173_c0"        # F1-ATPase (positive = synthesis, uses PMF)
NA_H_ANTI    = "rxn05209_c0"        # Na+/H+ antiporter
FLAGELLA_RXN = "flagella_Na_import_c0"

# Dissipation sweep: 0 to 100 in steps of 4
SWEEP_VALUES = list(range(0, 105, 4))

# Organic acid / overflow exchange reactions to track
OVERFLOW_RXNS = {
    "EX_cpd00029_e0": "acetate",
    "EX_cpd00020_e0": "pyruvate",
    "EX_cpd00036_e0": "succinate",
    "EX_cpd00047_e0": "formate",
    "EX_cpd00159_e0": "lactate",
    "EX_cpd00106_e0": "fumarate",
    "EX_cpd00013_e0": "NH3_secretion",   # nitrogen released from amino acids
}

# Substrate panel: (name, met_id, n_c, uses_Na_symporter, symporter_rxn_id)
# Two non-Na+ controls (glucose via H+ symport, glycerol via facilitated diffusion)
# Five Na+-symporter amino acids
# Acetate uses Na+ symport but is also an overflow metabolite -- included for
# completeness but acetate secretion obviously cannot be observed on acetate.
SUBSTRATES = [
    {"name": "glucose",   "met_id": "cpd00027", "n_c": 6, "na_symporter": False, "symporter_rxn": None},
    {"name": "glycerol",  "met_id": "cpd00100", "n_c": 3, "na_symporter": False, "symporter_rxn": None},
    {"name": "glutamate", "met_id": "cpd00023", "n_c": 5, "na_symporter": True,  "symporter_rxn": "rxn05298_c0"},
    {"name": "aspartate", "met_id": "cpd00041", "n_c": 4, "na_symporter": True,  "symporter_rxn": "rxn34493_c0"},
    {"name": "alanine",   "met_id": "cpd00035", "n_c": 3, "na_symporter": True,  "symporter_rxn": "rxn05215_c0"},
    {"name": "glycine",   "met_id": "cpd00033", "n_c": 2, "na_symporter": True,  "symporter_rxn": "rxn08661_c0"},
    {"name": "lysine",    "met_id": "cpd00039", "n_c": 6, "na_symporter": True,  "symporter_rxn": "rxn08854_c0"},
]


def add_flagella_reaction(model: cobra.Model) -> cobra.Reaction:
    """Add the flagellar Na+ import reaction to the model (in place)."""
    rxn = cobra.Reaction(FLAGELLA_RXN)
    rxn.name = "Flagella Na+ import (SMF dissipation)"
    rxn.add_metabolites({
        model.metabolites.cpd00971_e0: -1,
        model.metabolites.cpd00971_c0:  1,
    })
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    model.add_reactions([rxn])
    return rxn


def run_pfba_safe(model: cobra.Model) -> tuple[cobra.Solution | None, bool]:
    """Run pFBA, return (solution, feasible). Returns (None, False) on failure."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = cobra.flux_analysis.pfba(model)
        if sol.status == "optimal":
            return sol, True
        return None, False
    except Exception:
        return None, False


def extract_result_row(
    sol: cobra.Solution,
    substrate_name: str,
    dissipation_lb: float,
    symporter_rxn: str | None,
) -> dict:
    """Pull the metrics we care about from a pFBA solution."""
    row = {
        "substrate":       substrate_name,
        "dissipation_lb":  dissipation_lb,
        "feasible":        True,
        "growth_rate":     sol.fluxes[BIOMASS_RXN],
        "NaNQR_flux":      sol.fluxes[NANQR_RXN],
        "Na_EX_flux":      sol.fluxes[NA_EX_RXN],   # negative = importing from medium
        "ATPsyn_flux":     sol.fluxes[ATPSYN_RXN],
        "NaH_anti_flux":   sol.fluxes[NA_H_ANTI],
        "flagella_flux":   sol.fluxes[FLAGELLA_RXN],
        "symporter_flux":  sol.fluxes[symporter_rxn] if symporter_rxn else 0.0,
    }
    # Add overflow fluxes
    for rxn_id, name in OVERFLOW_RXNS.items():
        if rxn_id in sol.fluxes.index:
            row[name] = sol.fluxes[rxn_id]
        else:
            row[name] = 0.0
    return row


def infeasible_row(substrate_name: str, dissipation_lb: float) -> dict:
    row = {
        "substrate":       substrate_name,
        "dissipation_lb":  dissipation_lb,
        "feasible":        False,
        "growth_rate":     float("nan"),
        "NaNQR_flux":      float("nan"),
        "Na_EX_flux":      float("nan"),
        "ATPsyn_flux":     float("nan"),
        "NaH_anti_flux":   float("nan"),
        "flagella_flux":   float("nan"),
        "symporter_flux":  float("nan"),
    }
    for name in OVERFLOW_RXNS.values():
        row[name] = float("nan")
    return row


# ---------------------------------------------------------------------------
# Load model and media
# ---------------------------------------------------------------------------
print("Loading model...")
model_base = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
add_flagella_reaction(model_base)

media_defs = MEDIA
minimal_media = media_defs["minimal"]

# ---------------------------------------------------------------------------
# Baseline sanity check (zero dissipation on each substrate)
# ---------------------------------------------------------------------------
print("\n=== Baseline sanity check (no SMF dissipation) ===")
print(f"{'Substrate':<12} {'Growth':>8} {'NaNQR':>8} {'Na_EX':>8} {'Symporter':>10} {'Status'}")
print("-" * 60)

for sub in SUBSTRATES:
    with model_base:
        ex_id = f"EX_{sub['met_id']}_e0"
        media = minimal_media.copy()
        media[ex_id] = TOTAL_UPTAKE / sub["n_c"]
        model_base.medium = media
        model_base.reactions.get_by_id(FLAGELLA_RXN).lower_bound = 0.0

        sol, ok = run_pfba_safe(model_base)
        if ok:
            sym_flux = sol.fluxes[sub["symporter_rxn"]] if sub["symporter_rxn"] else 0.0
            print(
                f"{sub['name']:<12} {sol.fluxes[BIOMASS_RXN]:>8.4f}"
                f" {sol.fluxes[NANQR_RXN]:>8.3f}"
                f" {sol.fluxes[NA_EX_RXN]:>8.3f}"
                f" {sym_flux:>10.3f}  ok"
            )
        else:
            print(f"{sub['name']:<12}  INFEASIBLE at baseline!")

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
print("\n=== Running 2D sweep (substrate × dissipation level) ===")
all_results = []
# Track where each substrate first becomes infeasible
infeasibility_thresholds = {}

for sub in SUBSTRATES:
    ex_id = f"EX_{sub['met_id']}_e0"
    substrate_name = sub["name"]
    print(f"\n  {substrate_name} (n_c={sub['n_c']}, Na+-symporter={sub['na_symporter']})")

    # Baseline NaNQR flux (for context)
    with model_base:
        media = minimal_media.copy()
        media[ex_id] = TOTAL_UPTAKE / sub["n_c"]
        model_base.medium = media
        model_base.reactions.get_by_id(FLAGELLA_RXN).lower_bound = 0.0
        sol0, _ = run_pfba_safe(model_base)
        baseline_nanqr = sol0.fluxes[NANQR_RXN] if sol0 else float("nan")
        baseline_growth = sol0.fluxes[BIOMASS_RXN] if sol0 else float("nan")

    print(f"    baseline growth={baseline_growth:.4f}, NaNQR={baseline_nanqr:.2f}")

    hit_infeasible = False
    for lb in SWEEP_VALUES:
        with model_base:
            media = minimal_media.copy()
            media[ex_id] = TOTAL_UPTAKE / sub["n_c"]
            model_base.medium = media
            model_base.reactions.get_by_id(FLAGELLA_RXN).lower_bound = lb

            sol, ok = run_pfba_safe(model_base)

            if ok:
                row = extract_result_row(sol, substrate_name, lb, sub["symporter_rxn"])

                # Check Na+ trivially-free condition:
                # If EX_cpd00971_e0 flux < -1 (importing > 1 mmol/gDW/hr from medium),
                # the model is partly satisfying flagella demand from medium Na+ rather
                # than purely through NaNQR cycling. Flag this.
                if sol.fluxes[NA_EX_RXN] < -1.0:
                    row["na_medium_import_flag"] = True
                    print(
                        f"    lb={lb:4d}: WARNING Na+ medium import={sol.fluxes[NA_EX_RXN]:.1f}"
                        f" (cell drawing Na+ from medium)"
                    )
                else:
                    row["na_medium_import_flag"] = False

                all_results.append(row)

                # Save flux distribution if interesting:
                # (a) any overflow > 0.05, or (b) growth < 80% of baseline
                total_overflow = sum(
                    sol.fluxes.get(rxn_id, 0)
                    for rxn_id in OVERFLOW_RXNS
                    if rxn_id != "EX_cpd00013_e0"   # NH3 is expected on amino acids
                )
                growth_fraction = sol.fluxes[BIOMASS_RXN] / baseline_growth if baseline_growth > 0 else 1.0
                if total_overflow > 0.05 or growth_fraction < 0.80:
                    fname = FLUX_DIR / f"fluxes_{substrate_name}_lb{lb:04d}.json"
                    sol.fluxes.to_json(fname)

            else:
                if not hit_infeasible:
                    infeasibility_thresholds[substrate_name] = lb
                    print(f"    lb={lb:4d}: INFEASIBLE (threshold={lb})")
                    hit_infeasible = True
                all_results.append(infeasible_row(substrate_name, lb))

# ---------------------------------------------------------------------------
# Save main results
# ---------------------------------------------------------------------------
results_df = pd.DataFrame(all_results)
results_df.to_csv(OUT_PATH / "smf_sweep_results.csv", index=False)
print(f"\nResults saved to {OUT_PATH / 'smf_sweep_results.csv'}")

# ---------------------------------------------------------------------------
# Print infeasibility summary
# ---------------------------------------------------------------------------
print("\n=== Infeasibility thresholds (min dissipation level causing INFEASIBLE) ===")
for sub in SUBSTRATES:
    name = sub["name"]
    thresh = infeasibility_thresholds.get(name, ">100 (never infeasible in sweep)")
    print(f"  {name:<12}: {thresh}")

# ---------------------------------------------------------------------------
# Quick Na+ balance summary to check trivially-free scenario
# ---------------------------------------------------------------------------
print("\n=== Na+ medium import flag summary ===")
flagged = results_df[results_df.get("na_medium_import_flag", False) == True]
if len(flagged) == 0:
    print("  Na+ medium import never exceeded 1 mmol/gDW/hr.")
    print("  Na+ is cycling internally through NaNQR -- experiment is meaningful.")
else:
    print(f"  {len(flagged)} conditions show Na+ medium import > 1 mmol/gDW/hr.")
    print("  Check these conditions for interpretation.")
    print(flagged[["substrate", "dissipation_lb", "Na_EX_flux"]].to_string(index=False))

print("\nDone.")

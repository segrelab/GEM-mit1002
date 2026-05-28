#!/usr/bin/env python3
"""
Run diel FBA simulations for Alteromonas MIT1002 growth on Pro-released exometabolites.

For each time interval in the diel cycle and each value of f (the interception
fraction × Pro:Alt cell ratio), converts empirical per-Pro-cell release rates
to FBA exchange bounds and runs pFBA.

Ten of the twelve measured metabolites lack extracellular forms in iHS4156.
This script adds exchange + uptake transport reactions for those ten in-memory
only — model.xml is never modified. All new transport reactions are
irreversible (uptake-only): cpd_e0 --> cpd_c0.

Flux bound formula:
    bound [mmol/gDW/hr] = f × rate [nmol/Pro-cell/hr] × 1e-6 / alt_dw [g]

Alteromonas dry weight: 250 fg = 2.5e-13 g (placeholder; see Pedler et al.
2014 PNAS for MIT1002; flag for literature review).
"""

import warnings
from pathlib import Path

import cobra
import cobra.flux_analysis
import numpy as np
import pandas as pd

# ── Parameters ─────────────────────────────────────────────────────────────────

ALT_DW_G = 2.5e-13  # 250 fg, Alteromonas dry weight per cell

F_VALUES = [10.0]  # interception fraction × Pro:Alt cell ratio

# ── NGAM handling ───────────────────────────────────────────────────────────────
# The model's ATP maintenance reaction (rxn00062_c0, ATP hydrolysis) has lb=6.86
# mmol/gDW/hr — calibrated from lab batch-culture data. This is ~20× larger than
# the total carbon flux available from Pro exometabolites at f=10 (peak interval),
# making every simulation infeasible.
#
# For in-situ slow-growth conditions (μ ~ 0.001–0.01 h⁻¹), maintenance costs
# are far lower than lab-calibrated values. Setting NGAM_LB=0 allows the model
# to operate at any substrate level; pFBA will still minimise total flux, which
# implicitly penalises gratuitous ATP hydrolysis.
#
# FLAG FOR HELEN: This assumption needs explicit justification in the manuscript.
# Options: (1) keep 0 and cite the oligotrophic-condition precedent,
# (2) use a literature-derived low-maintenance value (~0.5 mmol/gDW/hr).
NGAM_RXN_ID = "rxn00062_c0"
NGAM_LB = 0.0  # set to 0 for in-situ slow-growth; original lab value = 6.86

# Per-solve timeout for GLPK (seconds). Some intervals are LP-degenerate at
# certain f values and cause the simplex method to hang; this caps wall time.
SOLVER_TIMEOUT_S = 30

# ── Paths ───────────────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).parent
MODEL_FILE = SCRIPT_DIR / "../../model.xml"
RATES_FILE = SCRIPT_DIR / "results/ProDiel_per_pro_cell_rates.csv"
MAP_FILE = SCRIPT_DIR / "metabolite_id_map.csv"
OUT_DIR = SCRIPT_DIR / "results"
OUT_DIR.mkdir(exist_ok=True)

# ── Basal inorganic medium ──────────────────────────────────────────────────────

BASAL_MEDIUM = {
    "EX_cpd00007_e0": 20,  # O2
    "EX_cpd00067_e0": 1000,  # H+  # FIXME: Should I remove this?
    # "EX_cpd00013_e0": 1000,  # NH3  (Remove to see if Amac can recycle N from Pro)
    "EX_cpd00058_e0": 1000,  # Cu2+
    "EX_cpd00971_e0": 1000,  # Na+
    "EX_cpd00063_e0": 1000,  # Ca2+
    "EX_cpd00048_e0": 1000,  # Sulfate
    "EX_cpd10516_e0": 1000,  # fe3
    "EX_cpd00254_e0": 1000,  # Mg
    "EX_cpd00009_e0": 1000,  # Phosphate
    "EX_cpd00205_e0": 1000,  # K+
    "EX_cpd00099_e0": 1000,  # Cl-
    "EX_cpd00030_e0": 1000,  # Mn2+
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00034_e0": 1000,  # Zn2+
    "EX_cpd00149_e0": 1000,  # Co2+
}

# ── Metabolites needing in-memory transport reactions ───────────────────────────
# These have cytosolic forms in iHS4156 but no extracellular metabolite,
# exchange reaction, or transport reaction.

MISSING_TRANSPORTS = {
    "cpd02711": "2-Keto-3-deoxy-6-phosphogluconate",
    "cpd00169": "3-Phosphoglycerate",
    "cpd00024": "2-Oxoglutarate (alpha-ketoglutarate)",
    "cpd00018": "AMP",
    "cpd00395": "Cysteate",
    "cpd01311": "Desthiobiotin",
    "cpd00111": "Oxidized glutathione (GSSG)",
    "cpd00477": "N-Acetylglutamate",
    "cpd00061": "Phosphoenolpyruvate",
    "cpd00091": "UMP",
}


def main() -> None:
    print("=" * 60)
    print("Diel FBA — Alteromonas MIT1002 on Pro exometabolites")
    print("=" * 60)

    print(f"\nAlteromonas dry weight: {ALT_DW_G:.2e} g ({ALT_DW_G * 1e15:.0f} fg)")
    print(f"f values: {F_VALUES}")
    print(
        f"NGAM lb: {NGAM_LB} mmol/gDW/hr (original lab value = 6.86; set to 0 for in-situ)"
    )

    # Load model
    print("\nLoading model (iHS4156)...")
    model = cobra.io.read_sbml_model(str(MODEL_FILE))
    print(f"  {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

    # Set ATP Maintenance to 0
    ngam_rxn = model.reactions.get_by_id(NGAM_RXN_ID)
    original_ngam_lb = ngam_rxn.lower_bound
    ngam_rxn.lower_bound = NGAM_LB
    print(f"\nNGAM ({NGAM_RXN_ID}): lb {original_ngam_lb} → {NGAM_LB}")

    # Configure GLPK: enable presolve and per-solve timeout.
    # Without presolve, GLPK simplex hangs indefinitely on some intervals due to
    # warm-start state from a previous LP interacting badly with very small bounds
    # (rates × f can produce bounds ~1e-7 mmol/gDW/hr). Timeout caps any runaway.
    model.solver.configuration.presolve = True
    model.solver.configuration.timeout = SOLVER_TIMEOUT_S
    print(f"GLPK: presolve=True, per-solve timeout={SOLVER_TIMEOUT_S}s")

    # Add transport and exchange reactions for all of the ProDiel metabolites
    print(
        "\nAdding transport reactions for metabolites absent from extracellular space..."
    )
    added = add_transport_reactions(model)
    print(
        f"  Model after additions: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites"
    )

    # Load rates and map
    print("\nLoading rates and metabolite→cpd map...")
    rates, met_map = load_inputs(RATES_FILE, MAP_FILE)
    n_intervals = rates.groupby(["interval_start_h", "interval_end_h"]).ngroups
    print(
        f"  {len(rates)} rows, {rates['metabolite'].nunique()} metabolites, {n_intervals} intervals"
    )

    # Run simulations
    print(
        f"\nRunning FBA: {n_intervals} intervals × {len(F_VALUES)} f values = {n_intervals * len(F_VALUES)} solves..."
    )
    growth_df, flux_df, flux_jsons, binding_df = run_simulations(
        model, rates, met_map, F_VALUES
    )

    # Save outputs
    growth_file = OUT_DIR / "growth_rates.csv"
    flux_file = OUT_DIR / "fluxes_long.csv"
    binding_file = OUT_DIR / "binding_constraints.csv"
    for key, flux_dict in flux_jsons.items():
        flux_json_file = OUT_DIR / f"fluxes_{key}.json"
        pd.Series(flux_dict).to_json(flux_json_file)

    growth_df.to_csv(growth_file, index=False)
    flux_df.to_csv(flux_file, index=False)
    binding_df.to_csv(binding_file, index=False)

    print("\nOutputs saved:")
    print(f"  {growth_file}  ({len(growth_df)} rows)")
    print(f"  {flux_file}  ({len(flux_df):,} rows)")
    print(f"  {binding_file}  ({len(binding_df)} rows)")

    # Summary table
    print("\nGrowth rate summary by f (across all intervals):")
    summary = (
        growth_df[growth_df["status"] == "optimal"]
        .groupby("f")["growth_rate"]
        .agg(mean="mean", min="min", max="max", n_feasible="count")
    )
    print(summary.to_string())
    print(f"\nTotal intervals: {n_intervals}  (infeasible at f=0 expected)")


# Helper functions
def add_transport_reactions(model: cobra.Model) -> list[str]:
    """Add exchange + uptake transport reactions for metabolites absent from
    the extracellular space. Modifies the model in-place; model.xml is unchanged.

    For each cpd in MISSING_TRANSPORTS:
      - Creates cpd_e0 metabolite (formula/charge copied from cytosolic form)
      - Adds exchange reaction EX_cpd_e0: cpd_e0 --> (lb=0, ub=1000)
        model.medium will set lb to a negative value to allow uptake
      - Adds uptake transport TP_cpd_e0: cpd_e0 --> cpd_c0 (irreversible, lb=0)

    Returns list of cpd IDs successfully added.
    """
    added = []
    for cpd_id, cpd_name in MISSING_TRANSPORTS.items():
        cyto_id = f"{cpd_id}_c0"
        ext_id = f"{cpd_id}_e0"
        ex_rxn_id = f"EX_{cpd_id}_e0"
        tp_rxn_id = f"TP_{cpd_id}_e0"

        if cyto_id not in model.metabolites:
            print(f"  WARNING: {cyto_id} not in model — skipping {cpd_name}")
            continue

        cyto_met = model.metabolites.get_by_id(cyto_id)

        ext_met = cobra.Metabolite(
            ext_id,
            name=f"{cpd_name} [e0]",
            compartment="e0",
            formula=cyto_met.formula,
            charge=cyto_met.charge,
        )

        # Exchange: cpd_e0 --> (stoich -1 on ext_met = consumed in forward = secretion)
        # Negative flux (uptake) is enabled when model.medium sets lb = -bound.
        ex_rxn = cobra.Reaction(
            ex_rxn_id,
            name=f"Exchange for {cpd_name} [e0]",
            lower_bound=0,
            upper_bound=1000,
        )
        ex_rxn.add_metabolites({ext_met: -1})

        # Uptake transport: cpd_e0 --> cpd_c0 (irreversible)
        tp_rxn = cobra.Reaction(
            tp_rxn_id,
            name=f"Uptake of {cpd_name}",
            lower_bound=0,
            upper_bound=1000,
        )
        tp_rxn.add_metabolites({ext_met: -1, cyto_met: 1})

        model.add_reactions([ex_rxn, tp_rxn])
        added.append(cpd_id)
        print(f"  Added: EX_{cpd_id}_e0  +  TP_{cpd_id}_e0  ({cpd_name})")

    return added


# ── Data loading ────────────────────────────────────────────────────────────────


def load_inputs(
    rates_file: Path, map_file: Path
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load per-Pro-cell release rates and metabolite name → cpd ID mapping."""
    rates = pd.read_csv(rates_file)
    met_map_df = pd.read_csv(map_file)
    met_map = dict(zip(met_map_df["name"], met_map_df["id"]))
    return rates, met_map


# ── Flux bound conversion ───────────────────────────────────────────────────────


def rate_to_bound(rate_nmol_per_hr: float, f: float) -> float:
    """Convert per-Pro-cell release rate to Alteromonas FBA exchange bound.

    bound [mmol gDW⁻¹ hr⁻¹] = f × rate [nmol cell⁻¹ hr⁻¹] × 1e-6 / alt_dw [g]
    """
    return f * rate_nmol_per_hr * 1e-6 / ALT_DW_G


# ── Simulation loop ─────────────────────────────────────────────────────────────


def run_simulations(
    model: cobra.Model,
    rates: pd.DataFrame,
    met_map: dict[str, str],
    f_values: list[float],
) -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    """Run FBA for every (interval, f) combination.

    Returns three tidy DataFrames:
      growth_rates:        interval_start_h, interval_end_h, f, growth_rate, status
      fluxes_long:         interval_start_h, interval_end_h, f, reaction_id, flux
      binding_constraints: interval_start_h, interval_end_h, f,
                           exchange_reaction, bound_value, flux_at_optimum, is_binding
    """
    intervals = rates.groupby(["interval_start_h", "interval_end_h"])
    total = intervals.ngroups * len(f_values)
    n = 0

    growth_rows: list[dict] = []
    flux_rows: list[dict] = []
    binding_rows: list[dict] = []

    for (t_start, t_end), interval_df in intervals:
        # Map metabolite names to exchange reaction IDs and their rates
        pro_rates: dict[str, float] = {}  # ex_rxn_id -> rate [nmol/Pro-cell/hr]
        for _, row in interval_df.iterrows():
            rate = row["per_pro_cell_rate_nmol_per_hr"]
            if rate <= 0:
                continue
            cpd_id = met_map.get(row["metabolite"])
            if cpd_id is None:
                continue
            ex_rxn_id = f"EX_{cpd_id}_e0"
            if ex_rxn_id not in model.reactions:
                continue
            pro_rates[ex_rxn_id] = rate

        for f in f_values:
            n += 1

            # Build medium: basal inorganics + Pro-released metabolites
            medium = BASAL_MEDIUM.copy()
            active_bounds: dict[str, float] = {}
            for ex_rxn_id, rate in pro_rates.items():
                bound = rate_to_bound(rate, f)
                if bound > 0:
                    medium[ex_rxn_id] = bound
                    active_bounds[ex_rxn_id] = bound

            # Run FBA (plain optimize; pFBA hangs on GLPK for some intervals)
            with model:
                model.medium = medium
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        solution = model.optimize()
                    if solution.status == "optimal":
                        growth_rate = float(solution.fluxes["bio1_biomass"])
                        status = "optimal"
                        fluxes = solution.fluxes
                    else:
                        growth_rate = float("nan")
                        status = solution.status
                        fluxes = None
                except Exception:
                    growth_rate = float("nan")
                    status = "infeasible"
                    fluxes = None

            mu_str = f"{growth_rate:.5f}" if not np.isnan(growth_rate) else "infeasible"
            print(
                f"  [{n:3d}/{total}]  t={t_start:2g}-{t_end:2g}h  f={f:4g}  μ={mu_str}"
            )

            growth_rows.append(
                {
                    "interval_start_h": t_start,
                    "interval_end_h": t_end,
                    "f": f,
                    "growth_rate": growth_rate,
                    "status": status,
                }
            )

            if fluxes is not None:
                for rxn_id, flux in fluxes.items():
                    flux_rows.append(
                        {
                            "interval_start_h": t_start,
                            "interval_end_h": t_end,
                            "f": f,
                            "reaction_id": rxn_id,
                            "flux": flux,
                        }
                    )

                # Binding constraints: which uptake bounds are active at the optimum?
                # Uptake exchange runs at negative flux; |flux| ≈ bound when binding.
                for ex_rxn_id, bound in active_bounds.items():
                    if ex_rxn_id not in fluxes.index:
                        continue
                    actual_flux = float(fluxes[ex_rxn_id])
                    binding_rows.append(
                        {
                            "interval_start_h": t_start,
                            "interval_end_h": t_end,
                            "f": f,
                            "exchange_reaction": ex_rxn_id,
                            "bound_value": bound,
                            "flux_at_optimum": actual_flux,
                            "is_binding": abs(actual_flux) >= 0.99 * bound,
                        }
                    )

    # Convery flux rows into a pandas dataframe
    flux_rows = pd.DataFrame(flux_rows)

    # Make a dictionary of fluxes in JSON-serializable format for quick lookup in plotting
    flux_jsons = {}
    for (t_start, t_end, f), group in flux_rows.groupby(
        ["interval_start_h", "interval_end_h", "f"]
    ):
        key = f"{t_start}_{t_end}_{f}"
        flux_jsons[key] = dict(zip(group["reaction_id"], group["flux"]))

    return (
        pd.DataFrame(growth_rows),
        flux_rows,
        flux_jsons,
        pd.DataFrame(binding_rows),
    )


if __name__ == "__main__":
    main()

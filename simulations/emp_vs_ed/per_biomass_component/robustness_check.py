"""Robustness check on the ED usage seen in the per-biomass-component maximization.

For each amino acid that showed ED flux when its production was maximized, we ask
three questions:

  1. Is ED *required* to achieve maximal production, or is it just one of several
     equally-optimal solutions that pFBA happens to pick?  -> FVA on the ED
     reaction at (near-)optimal production. If the FVA minimum is > 0, ED is
     obligately required; if it is 0, ED is an interchangeable alternative.

  2. Does the ED usage survive loopless FBA (i.e. is it not an artifact of a
     thermodynamically-infeasible cycle)? -> loopless FVA on the ED reaction.

  3. What is ED providing? -> at a representative pFBA solution, report the
     ED-coupled NADPH source (G6PDH) and the competing oxidative-PPP branch, so
     we can tell a redox story from a carbon (pyruvate) story.

Setup mirrors per_biomass_component.ipynb: glucose minimal medium, sink reactions
(lb=0) added for every metabolite, objective = maximize the sink of the target.
"""

from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, loopless_solution, pfba

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
OUT_PATH = FILE_PATH / "robustness_results"
OUT_PATH.mkdir(exist_ok=True)

# Marker reactions (consistent with the other emp_vs_ed scripts)
ED_RXN = "rxn01477_c0"       # KDPG aldolase (eda) -- ED marker
EMP_RXN = "rxn00558_c0"      # PFK -- EMP marker
G6PDH_RXN = "rxn00604_c0"    # glucose-6-P dehydrogenase (NADPH; shared ED/oxPPP entry)
OXPPP_RXN = "rxn01115_c0"    # 6-phosphogluconate dehydrogenase (oxPPP branch, NADPH + CO2)
BIOMASS_RXN = "bio1_biomass"

# Amino acids to test (name -> cytosolic metabolite id)
AMINO_ACIDS = {
    "L-Cysteine": "cpd00084_c0",
    "L-Leucine": "cpd00107_c0",
    "L-Glutamate": "cpd00023_c0",
    "L-Glutamine": "cpd00053_c0",
    "L-Alanine": "cpd00035_c0",
    "L-Valine": "cpd00156_c0",
    "L-Proline": "cpd00129_c0",
}

# Medium from per_biomass_component.ipynb (glucose minimal)
MEDIUM = {
    "EX_cpd00027_e0": 10,    # Glucose
    "EX_cpd00007_e0": 20,    # Oxygen
    "EX_cpd00013_e0": 1000,  # Ammonia
    "EX_cpd00011_e0": 1000,  # CO2
    "EX_cpd00067_e0": 1000,  # H+
    "EX_cpd00009_e0": 1000,  # Phosphate
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00063_e0": 1000,  # Ca2+
    "EX_cpd00099_e0": 1000,  # Cl-
    "EX_cpd00149_e0": 1000,  # Co2+
    "EX_cpd00058_e0": 1000,  # Cu2+
    "EX_cpd00254_e0": 1000,  # Mg2+
    "EX_cpd00205_e0": 1000,  # K+
    "EX_cpd00971_e0": 1000,  # Na+
    "EX_cpd00048_e0": 1000,  # Sulfate
    "EX_cpd00034_e0": 1000,  # Zn2+
    "EX_cpd10516_e0": 1000,  # Fe+3
    "EX_cpd00030_e0": 1000,  # Mn2+
}

FRACTION = 0.999  # near-optimal production for FVA
TOL = 1e-6        # treat |flux| below this as zero


def build_model():
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
    model.medium = MEDIUM
    # Add sink reactions (lb=0) for every metabolite, as in the notebook
    existing = {r.id for r in model.reactions}
    for met in model.metabolites:
        if "SK_" + met.id not in existing:
            model.add_boundary(met, type="sink", lb=0)
    return model


def analyze(model, name, met_id):
    sink_id = "SK_" + met_id
    rec = {"amino_acid": name}
    with model:
        model.objective = {model.reactions.get_by_id(sink_id): 1}

        # Max production and a representative parsimonious solution
        fba = model.optimize()
        max_prod = fba.objective_value
        rec["max_production"] = max_prod
        if fba.status != "optimal" or max_prod < TOL:
            rec["feasible"] = False
            return rec
        rec["feasible"] = True

        sol = pfba(model)
        rec["ed_flux_pfba"] = sol.fluxes[ED_RXN]
        rec["emp_flux_pfba"] = sol.fluxes[EMP_RXN]
        rec["g6pdh_flux_pfba"] = sol.fluxes[G6PDH_RXN]
        rec["oxppp_flux_pfba"] = sol.fluxes[OXPPP_RXN]

        # Loopless projection of that solution
        ll = loopless_solution(model)
        rec["ed_flux_loopless"] = ll.fluxes[ED_RXN]

        # Standard FVA: is ED required at (near-)optimal production?
        fva = flux_variability_analysis(
            model, reaction_list=[ED_RXN, EMP_RXN], fraction_of_optimum=FRACTION
        )
        rec["ed_min"] = fva.loc[ED_RXN, "minimum"]
        rec["ed_max"] = fva.loc[ED_RXN, "maximum"]
        rec["emp_min"] = fva.loc[EMP_RXN, "minimum"]
        rec["emp_max"] = fva.loc[EMP_RXN, "maximum"]

        # Loopless FVA on ED: required even with cycles forbidden?
        try:
            fva_ll = flux_variability_analysis(
                model, reaction_list=[ED_RXN], fraction_of_optimum=FRACTION,
                loopless=True,
            )
            rec["ed_min_loopless"] = fva_ll.loc[ED_RXN, "minimum"]
            rec["ed_max_loopless"] = fva_ll.loc[ED_RXN, "maximum"]
        except Exception as e:
            rec["ed_min_loopless"] = None
            rec["ed_max_loopless"] = None
            print(f"  loopless FVA failed for {name}: {e}")

    # Verdict
    rec["ed_required"] = bool(rec.get("ed_min", 0) > TOL)
    rec["ed_required_loopless"] = bool((rec.get("ed_min_loopless") or 0) > TOL)
    return rec


def main():
    model = build_model()
    records = [analyze(model, name, mid) for name, mid in AMINO_ACIDS.items()]
    for r in records:
        print(
            f"{r['amino_acid']:14s} feasible={r.get('feasible')} "
            f"ED_pfba={r.get('ed_flux_pfba', float('nan')):.3f} "
            f"ED_min={r.get('ed_min', float('nan')):.3f} "
            f"ED_min_loopless={r.get('ed_min_loopless')} "
            f"required={r.get('ed_required')} required_loopless={r.get('ed_required_loopless')}"
        )
    df = pd.DataFrame.from_records(records)
    df.to_csv(OUT_PATH / "aa_ed_robustness.csv", index=False)
    print("\nSaved", OUT_PATH / "aa_ed_robustness.csv")


if __name__ == "__main__":
    main()

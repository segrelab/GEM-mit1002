"""Evaluate the known growth phenotypes against the model.

Both the CI regression report (``scripts/generate_growth_report.py``) and the
manuscript figure (``scripts/generate_phenotype_figure.py``) need the same
thing: for every row of ``data/known_growth_phenotypes.tsv``, build the
medium, run FBA, and decide whether the prediction agrees with the
experiment. Keeping that here means the table and the figure cannot drift
apart, and a fix lands in both at once.

Three things this does that the loop inlined in ``generate_growth_report.py``
does not:

1. Metabolites are added to the medium one at a time. The original guarded
   the whole row with ``if all(exchange present for every metabolite)``, so a
   row like ``Methionine, Pyruvate`` had *neither* metabolite added when only
   methionine lacked an exchange reaction. The result is indistinguishable
   from "the model cannot grow on pyruvate", which silently turns a
   reconstruction gap into an apparent false negative.

2. Rows the model cannot represent are reported as ``no_exchange`` instead of
   as a prediction. A metabolite with no exchange reaction cannot enter the
   cell, so FBA returns zero growth for a trivial reason. Scoring those as
   correct "No" predictions inflates specificity: as of this writing 19 of
   the 61 rows are in this state and 14 of them are experimental "No".

3. Solver output is validated. Biomass has a lower bound of 0, so a negative
   objective value is not a growth rate. Those rows are reported as
   ``invalid_solve`` rather than plotted as no-growth.
"""

import os

import numpy as np
import pandas as pd

from tools.media import MEDIA

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Default location of the curated phenotype table.
PHENOTYPE_TSV = os.path.join(PROJECT_ROOT, "data", "known_growth_phenotypes.tsv")

#: Total carbon made available, in mmol C / gDW / hr. Uptake of each carbon
#: source is set to this divided by its number of carbons, so every condition
#: offers the same amount of carbon regardless of substrate. 60 matches a
#: glucose uptake of 10 mmol / gDW / hr.
TOTAL_C_UPTAKE = 60.0

#: Uptake bound for metabolites supplied as a nutrient other than carbon.
NON_CARBON_UPTAKE = 1000.0

#: Growth rates at or below this are called "no growth" (1 / hr).
GROWTH_THRESHOLD = 1e-3

#: Predictions below this multiple of the threshold are flagged as marginal:
#: they are the cells most likely to flip on an unrelated curation change.
NEAR_THRESHOLD_FACTOR = 10.0

#: Compounds supplied as nitrogen sources rather than as carbon sources. These
#: get a generous flat bound instead of the carbon-normalised one, so that a
#: nitrogen source which happens to contain carbon (urea) is not also dosed as
#: the condition's carbon source.
N_SOURCE_IDS = {
    "cpd00013": "Ammonium",
    "cpd00209": "Nitrate",
    "cpd00073": "Urea",
}

#: Compounds supplied as the carbon source in the nitrogen-source screen.
C_SOURCE_IDS = {
    "cpd00027": "Glucose",
    "cpd00020": "Pyruvate",
    "cpd00036": "Succinate",
}

CONCORDANT = ("true_positive", "true_negative")
DISCORDANT = ("false_positive", "false_negative")
#: Rows that carry no information about model quality either way.
UNSCORED = ("unsure", "no_exchange", "invalid_solve")


def exchange_id(met_id: str) -> str:
    """Exchange reaction id for a ModelSEED compound id."""
    return f"EX_{met_id}_e0"


def load_phenotypes(path: str | None = None) -> pd.DataFrame:
    """Read the phenotype table, splitting the comma-separated metabolite ids."""
    table = pd.read_csv(
        path or PHENOTYPE_TSV,
        sep="\t",
        converters={"met_id": lambda x: [m.strip() for m in x.split(",") if m.strip()]},
    )
    # Several rows have trailing whitespace in the display name ("Glucose ").
    for column in ("c_source", "minimal_media", "growth"):
        table[column] = table[column].astype(str).str.strip()
    return table


def _uptake_bound(model, met_id: str) -> float:
    """Carbon-normalised uptake bound for one metabolite."""
    if met_id in N_SOURCE_IDS:
        return NON_CARBON_UPTAKE
    try:
        metabolite = model.metabolites.get_by_id(met_id + "_c0")
    except KeyError:
        return NON_CARBON_UPTAKE
    n_carbon = metabolite.elements.get("C", 0)
    if n_carbon > 0:
        return TOTAL_C_UPTAKE / n_carbon
    return NON_CARBON_UPTAKE


def _solve(model) -> tuple[float, bool]:
    """Optimise and return ``(growth_rate, is_valid)``.

    A non-optimal status, a missing objective value, or a negative objective
    value all mean the number is not a growth rate. Biomass cannot carry
    negative flux, so a negative optimum indicates a solver or bounds problem
    rather than a phenotype.
    """
    solution = model.optimize()
    if solution.status != "optimal":
        return float("nan"), False
    value = solution.objective_value
    if value is None or not np.isfinite(value):
        return float("nan"), False
    if value < -1e-6:
        return float("nan"), False
    return max(float(value), 0.0), True


def _classify(experimental: str, predicted_growth: bool | None) -> str:
    if predicted_growth is None:
        return "no_exchange"
    if experimental not in ("Yes", "No"):
        return "unsure"
    observed_growth = experimental == "Yes"
    if observed_growth and predicted_growth:
        return "true_positive"
    if not observed_growth and not predicted_growth:
        return "true_negative"
    if predicted_growth:
        return "false_positive"
    return "false_negative"


def evaluate_phenotypes(
    model,
    phenotypes: pd.DataFrame | None = None,
    growth_threshold: float = GROWTH_THRESHOLD,
) -> pd.DataFrame:
    """Run every phenotype condition and classify the outcome.

    Returns a copy of the phenotype table with these columns added:

    ``missing_exchanges``
        Comma-separated compound ids that have no exchange reaction. A
        non-empty value means the condition could not be represented.
    ``evaluable``
        False when any required exchange reaction is absent, or when the
        solver returned something that is not a growth rate.
    ``fba_growth_rate``
        The growth rate, or NaN when the row is not evaluable.
    ``predicted``
        "Yes", "No", or None when not evaluable.
    ``category``
        One of true_positive, true_negative, false_positive, false_negative,
        unsure, no_exchange, invalid_solve.
    ``discordant``
        True for false positives and false negatives only.
    ``near_threshold``
        True when the predicted rate is positive but within
        ``NEAR_THRESHOLD_FACTOR`` of the threshold, so the call is marginal.
    """
    # Imported here rather than at module scope so that the plotting code can
    # import this module (for the category names and thresholds) without
    # needing cobra and GEM-utilities installed.
    from gem_utilities import media as media_utils

    if phenotypes is None:
        phenotypes = load_phenotypes()

    model_reactions = {reaction.id for reaction in model.reactions}
    records = []

    for _, row in phenotypes.iterrows():
        met_ids = list(row["met_id"])
        missing = [m for m in met_ids if exchange_id(m) not in model_reactions]
        present = [m for m in met_ids if exchange_id(m) in model_reactions]

        medium = MEDIA[row["minimal_media"]].copy()
        # Add every metabolite we can, even when a sibling is missing, so the
        # recorded growth rate reflects what the model actually did rather
        # than an empty medium.
        for met_id in present:
            medium[exchange_id(met_id)] = _uptake_bound(model, met_id)

        with model:
            model.medium = media_utils.clean_media(model, medium)
            rate, valid = _solve(model)

        if missing:
            category = "no_exchange"
            predicted = None
        elif not valid:
            category = "invalid_solve"
            predicted = None
        else:
            predicted_growth = rate > growth_threshold
            predicted = "Yes" if predicted_growth else "No"
            category = _classify(row["growth"], predicted_growth)

        records.append(
            {
                "missing_exchanges": ", ".join(missing),
                "evaluable": not missing and valid,
                "fba_growth_rate": rate if valid else float("nan"),
                "predicted": predicted,
                "category": category,
                "discordant": category in DISCORDANT,
                "near_threshold": bool(
                    valid
                    and 0 < rate <= NEAR_THRESHOLD_FACTOR * growth_threshold
                ),
            }
        )

    return pd.concat(
        [phenotypes.reset_index(drop=True), pd.DataFrame.from_records(records)],
        axis=1,
    )


def summarise(results: pd.DataFrame) -> dict:
    """Sensitivity and specificity over the rows that are actually scored.

    Reported separately on purpose. The dataset is roughly 80% growth-positive,
    so a single accuracy figure is dominated by true positives and says very
    little about whether the model can predict a *failure* to grow.
    """
    counts = results["category"].value_counts().to_dict()
    true_positive = counts.get("true_positive", 0)
    true_negative = counts.get("true_negative", 0)
    false_positive = counts.get("false_positive", 0)
    false_negative = counts.get("false_negative", 0)

    observed_growth = true_positive + false_negative
    observed_no_growth = true_negative + false_positive

    return {
        "n_conditions": len(results),
        "n_no_exchange": counts.get("no_exchange", 0),
        "n_invalid_solve": counts.get("invalid_solve", 0),
        "n_unsure": counts.get("unsure", 0),
        "n_scored": true_positive + true_negative + false_positive + false_negative,
        "true_positive": true_positive,
        "true_negative": true_negative,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "sensitivity": true_positive / observed_growth if observed_growth else np.nan,
        "specificity": (
            true_negative / observed_no_growth if observed_no_growth else np.nan
        ),
    }


def format_summary(summary: dict) -> str:
    """One-paragraph text version of :func:`summarise`, for logs and captions."""
    return (
        f"{summary['n_conditions']} conditions; "
        f"{summary['n_no_exchange']} not representable (no exchange reaction), "
        f"{summary['n_invalid_solve']} invalid solves, "
        f"{summary['n_unsure']} experimentally unsure; "
        f"{summary['n_scored']} scored. "
        f"Sensitivity {summary['sensitivity']:.2f} "
        f"({summary['true_positive']}/"
        f"{summary['true_positive'] + summary['false_negative']}), "
        f"specificity {summary['specificity']:.2f} "
        f"({summary['true_negative']}/"
        f"{summary['true_negative'] + summary['false_positive']})."
    )

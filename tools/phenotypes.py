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
UNSCORED = ("unsure", "no_exchange", "invalid_solve", "excluded")

#: Column naming a reason the row is not scored. Empty means "score this row".
EXCLUSION_COLUMN = "exclude_reason"

#: Closed vocabulary for :data:`EXCLUSION_COLUMN`. Keep in sync with the table
#: in ``data/README.md``; ``test_phenotype_data.py`` enforces it.
#:
#: The split that matters is between the first three and the last. The first
#: three are problems with the observation -- a better experiment could fix
#: them, and the row might come back. ``not_representable`` is a trusted
#: observation that a stoichiometric model cannot express even in principle,
#: which is a permanent property of the modelling formalism rather than a data
#: problem, and is worth reporting rather than quietly dropping.
EXCLUSION_REASONS = {
    "control_failed": (
        "The experiment's own control did not behave as required, so the "
        "condition cannot be interpreted."
    ),
    "id_uncertain": (
        "The compound could not be confidently mapped to a model metabolite."
    ),
    "conflicting_reports": (
        "The same condition was scored differently by different sources and "
        "the disagreement is unresolved."
    ),
    "not_representable": (
        "A trusted result that flux balance analysis cannot reproduce in "
        "principle -- regulation, inhibition, or a kinetic effect."
    ),
}

#: Baseline of mismatches that are known and accepted at the current state of
#: curation. See :func:`load_expected_mismatches`.
EXPECTED_MISMATCHES_TSV = os.path.join(
    PROJECT_ROOT, "test", "test_files", "expected_phenotype_mismatches.tsv"
)

MISMATCH_COLUMNS = ["minimal_media", "c_source", "category", "notes"]


def exchange_id(met_id: str) -> str:
    """Exchange reaction id for a ModelSEED compound id."""
    return f"EX_{met_id}_e0"


def condition_key(minimal_media: str, c_source: str) -> str:
    """Stable identifier for one row of the phenotype table.

    Medium plus displayed substrate is unique across the table and stays
    readable in a file people edit by hand, which matters more here than
    robustness -- ``test_phenotype_data.py`` fails if the pair stops being
    unique.
    """
    return f"{str(minimal_media).strip()} | {str(c_source).strip()}"


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
    # Tolerate the column being absent so older copies of the file still load.
    if EXCLUSION_COLUMN not in table.columns:
        table[EXCLUSION_COLUMN] = ""
    table[EXCLUSION_COLUMN] = (
        table[EXCLUSION_COLUMN].fillna("").astype(str).str.strip()
    )
    table["condition"] = [
        condition_key(m, c) for m, c in zip(table["minimal_media"], table["c_source"])
    ]
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
        unsure, no_exchange, invalid_solve, excluded. This is the single field
        downstream code should switch on.
    ``raw_category``
        The verdict the row would have had if it were not excluded. Identical
        to ``category`` for rows that are scored. Kept so that excluding a row
        never destroys the underlying comparison.
    ``discordant``
        True for false positives and false negatives that are actually scored.
        Always False for excluded rows, since an uninterpretable condition is
        not evidence that the model is wrong.
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
            raw_category = "no_exchange"
            predicted = None
        elif not valid:
            raw_category = "invalid_solve"
            predicted = None
        else:
            predicted_growth = rate > growth_threshold
            predicted = "Yes" if predicted_growth else "No"
            raw_category = _classify(row["growth"], predicted_growth)

        excluded = bool(row.get(EXCLUSION_COLUMN, ""))
        records.append(
            {
                "missing_exchanges": ", ".join(missing),
                "evaluable": not missing and valid,
                "fba_growth_rate": rate if valid else float("nan"),
                "predicted": predicted,
                "raw_category": raw_category,
                "category": "excluded" if excluded else raw_category,
                "discordant": (not excluded) and raw_category in DISCORDANT,
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
        "n_excluded": counts.get("excluded", 0),
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
        f"{summary['n_unsure']} experimentally unsure, "
        f"{summary['n_excluded']} excluded; "
        f"{summary['n_scored']} scored. "
        f"Sensitivity {summary['sensitivity']:.2f} "
        f"({summary['true_positive']}/"
        f"{summary['true_positive'] + summary['false_negative']}), "
        f"specificity {summary['specificity']:.2f} "
        f"({summary['true_negative']}/"
        f"{summary['true_negative'] + summary['false_positive']})."
    )


# --------------------------------------------------------------------------
# Expected-mismatch baseline
# --------------------------------------------------------------------------
#
# A model under curation always has mismatches; that is the normal state, not
# an error. A test that fails whenever any condition disagrees can only pass
# when the model is perfect, so it fails continuously and stops being read.
#
# Instead the currently-accepted mismatches are recorded in a committed file
# and the test compares against it. A *new* mismatch fails, and so does a
# listed mismatch that has started passing -- the second case matters just as
# much, because otherwise an accidental fix goes unnoticed and the baseline
# silently rots into a list of things that are no longer true.


def load_expected_mismatches(path: str | None = None) -> pd.DataFrame:
    """Read the accepted-mismatch baseline. A missing file reads as empty."""
    path = path or EXPECTED_MISMATCHES_TSV
    if not os.path.exists(path):
        return pd.DataFrame(columns=MISMATCH_COLUMNS + ["condition"])
    table = pd.read_csv(path, sep="\t").fillna("")
    for column in MISMATCH_COLUMNS:
        if column not in table.columns:
            table[column] = ""
    table["condition"] = [
        condition_key(m, c) for m, c in zip(table["minimal_media"], table["c_source"])
    ]
    return table


def write_expected_mismatches(
    results: pd.DataFrame, path: str | None = None, notes: str = ""
) -> pd.DataFrame:
    """Overwrite the baseline with the mismatches in ``results``.

    Deliberately not called from the test suite. Regenerating the baseline is
    a decision to accept the current state, so it belongs in a script someone
    runs on purpose and commits.
    """
    path = path or EXPECTED_MISMATCHES_TSV
    mismatches = results[results["discordant"]].copy()
    out = pd.DataFrame(
        {
            "minimal_media": mismatches["minimal_media"],
            "c_source": mismatches["c_source"],
            "category": mismatches["category"],
            "notes": notes,
        }
    ).sort_values(["minimal_media", "c_source"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    out.to_csv(path, sep="\t", index=False, lineterminator="\n")
    return out


def compare_to_baseline(
    results: pd.DataFrame, expected: pd.DataFrame | None = None
) -> dict:
    """Diff the current mismatches against the accepted baseline.

    Returns ``new`` (mismatches not in the baseline), ``resolved`` (baseline
    entries that now agree, or that no longer exist in the phenotype table),
    and ``changed`` (still mismatching, but as a different category -- a false
    positive that became a false negative is worth noticing).
    """
    if expected is None:
        expected = load_expected_mismatches()

    current = results[results["discordant"]]
    current_by_condition = dict(zip(current["condition"], current["category"]))
    expected_by_condition = dict(zip(expected["condition"], expected["category"]))
    known_conditions = set(results["condition"])

    new = sorted(set(current_by_condition) - set(expected_by_condition))
    resolved = sorted(set(expected_by_condition) - set(current_by_condition))
    changed = sorted(
        condition
        for condition in set(current_by_condition) & set(expected_by_condition)
        if current_by_condition[condition] != expected_by_condition[condition]
    )
    return {
        "new": [(c, current_by_condition[c]) for c in new],
        "resolved": [
            (c, expected_by_condition[c], c in known_conditions) for c in resolved
        ],
        "changed": [
            (c, expected_by_condition[c], current_by_condition[c]) for c in changed
        ],
    }

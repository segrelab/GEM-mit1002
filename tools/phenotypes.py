"""Evaluate the known growth phenotypes against the model.

Several places need the same thing: for every row of
``data/known_growth_phenotypes.tsv``, build the medium, run FBA, and decide
whether the prediction agrees with the experiment. Keeping it here means they
cannot drift apart, and a fix lands in all of them at once.

Current callers:

* ``test/test_growth.py`` -- the CI regression test against the accepted
  mismatch baseline.
* ``curation_process/run_tests_on_prs.py`` -- the same evaluation replayed
  across every merged PR, which is what figure 2B of the manuscript plots.

``scripts/generate_growth_report.py`` still carries its own inlined copy of
this loop and therefore still has all three bugs listed below. It only feeds a
human-read heatmap, not a number in the paper, but it should be migrated.

Three things this does that the loop inlined in ``generate_growth_report.py``
does not:

1. Metabolites are added to the medium one at a time. The original guarded
   the whole row with ``if all(exchange present for every metabolite)``, so a
   row like ``Methionine, Pyruvate`` had *neither* metabolite added when only
   methionine lacked an exchange reaction. The result is indistinguishable
   from "the model cannot grow on pyruvate", which silently turns a
   reconstruction gap into an apparent false negative.

2. An infeasible solve is read as no growth, not as a failure. See
   :func:`_solve`.

3. Solver output is validated. Biomass has a lower bound of 0, so a negative
   objective value is not a growth rate. Those rows are reported as
   ``invalid_solve`` rather than plotted as no-growth.

Missing exchange reactions
--------------------------

A condition whose compound has no exchange reaction in the model **is scored**,
as the no-growth prediction it is. An earlier version of this module held those
rows out as a separate ``no_exchange`` category on the grounds that FBA returns
no growth for a trivial reason.

That was wrong, for two independent reasons.

The first is curatorial. For most of these compounds the absence of a
transporter is a *finding*: the genome was searched, no candidate transporter
was found, and the model reflects that. "No genomic evidence for uptake,
therefore no growth" is a mechanistic prediction, and it is one of the few
kinds of prediction a stoichiometric model makes about *failure* to grow.
Discarding it removed most of the model's negative predictions and left
specificity resting on two conditions.

The second is arithmetic, and it is the reason this is scored on the model's
actual solution rather than by the shortcut of calling every missing exchange a
"No". Because point 1 above adds whatever metabolites it can, a multi-compound
condition can still grow on the compounds that are present. Both of these are
real rows:

* ``marine_broth_wo_yeast_and_peptone | Methionine, Pyruvate`` -- grows
  experimentally, has no methionine exchange, and the model grows on the
  pyruvate. Scoring it on the missing exchange would record a false negative
  for a condition the model gets right.
* ``marine_broth_wo_yeast_and_peptone | Cystine, Pyruvate`` -- does *not* grow
  experimentally, has no cystine exchange, and the model grows on the pyruvate
  anyway. Scoring it on the missing exchange would record a true negative for a
  condition the model gets wrong.

So the missing exchange is not itself the prediction; the solve is. Rows where
the no-growth call did rest on a missing uptake route are flagged in the
``no_uptake_route`` column, so a caption can report how many of the negative
predictions come from an absent transporter rather than from network structure.

The genuinely unscorable case still has a home: a trusted observation that FBA
cannot reproduce in principle goes in ``exclude_reason`` as
``not_representable``. That is a deliberate, documented, per-row judgement
rather than a side effect of the reconstruction's coverage.
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

#: Flux magnitude (mmol / gDW / hr) above which a reaction is treated as
#: implausibly large. Nothing in this organism should carry an order of
#: magnitude more flux than the carbon supply, so a reaction above this is
#: almost always a thermodynamically infeasible loop rather than biology.
#: Used by :func:`evaluate_phenotypes` when ``flux_limit`` is requested and by
#: ``curation_process/run_tests_on_prs.py``.
DEFAULT_FLUX_LIMIT = 100.0

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
#:
#: Note what is *not* here: a condition whose compound has no exchange reaction
#: in the model. That used to be held out as ``no_exchange``. It is now scored
#: as the no-growth prediction it is -- see the module docstring.
UNSCORED = ("unsure", "invalid_solve", "excluded")

#: Every value ``category`` can take. These partition the table: each row gets
#: exactly one, so the counts always sum to the number of phenotype rows.
#: ``test_phenotype_data.py`` enforces that, because a figure caption that says
#: "54 growth phenotypes" is only defensible if the categories add up.
CATEGORIES = CONCORDANT + DISCORDANT + UNSCORED

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

    An **infeasible** problem is a growth rate of zero, not an invalid result.
    The model carries a forced lower bound on ATP hydrolysis, so a medium that
    cannot supply maintenance energy has no feasible solution at all. That is
    not a solver failure; it is the model saying the organism cannot sustain
    itself on this medium, which is exactly the prediction the phenotype is
    being compared against. Reporting it as unscorable instead threw away real
    no-growth predictions -- most of them, in fact, since a medium whose only
    carbon source cannot be taken up is infeasible rather than zero-growth.

    Any *other* non-optimal status (unbounded, solver error), a missing
    objective value, or a negative one really is unusable. Biomass cannot carry
    negative flux, so a negative optimum indicates a solver or bounds problem
    rather than a phenotype.
    """
    solution = model.optimize()
    if solution.status == "infeasible":
        return 0.0, True
    if solution.status != "optimal":
        return float("nan"), False
    value = solution.objective_value
    if value is None or not np.isfinite(value):
        return float("nan"), False
    if value < -1e-6:
        return float("nan"), False
    return max(float(value), 0.0), True


def _large_flux_reactions(model, flux_limit: float) -> frozenset:
    """Reaction ids carrying more than ``flux_limit`` in a parsimonious solution.

    Plain FBA is useless for this question. Any thermodynamically infeasible
    loop in the network can carry unlimited flux without changing the objective,
    so the loop's magnitude in an FBA solution is whatever the solver happened
    to land on -- it can be enormous in one solve and zero in the next for an
    unchanged model. pFBA minimises total flux subject to the same optimum, so a
    reaction that still carries a huge flux is carrying it because the network
    forces it to, which is the thing worth reporting.

    An infeasible problem returns the empty set rather than raising. The model
    carries a non-zero maintenance requirement, so a medium that cannot supply
    that energy makes the LP infeasible; biologically that is just "no growth",
    and a condition with no growth has no fluxes to report.
    """
    # Imported here for the same reason as the media import below: this module
    # is imported by plotting code that has no cobra installed.
    import cobra

    try:
        solution = cobra.flux_analysis.pfba(model)
    except (cobra.exceptions.Infeasible, cobra.exceptions.OptimizationError):
        return frozenset()
    return frozenset(
        reaction_id
        for reaction_id, flux in solution.fluxes.items()
        if abs(flux) > flux_limit
    )


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
    flux_limit: float | None = None,
) -> pd.DataFrame:
    """Run every phenotype condition and classify the outcome.

    Pass ``flux_limit`` to also record, for every condition that grows, which
    reactions carry more flux than that. This costs a second (pFBA) solve per
    growing condition, so it is off by default: the test suite only needs the
    verdicts, while ``run_tests_on_prs.py`` needs both metrics and must get them
    from the same simulation to be comparable.

    Returns a copy of the phenotype table with these columns added:

    ``missing_exchanges``
        Comma-separated compound ids that have no exchange reaction. The
        condition is still scored; see the module docstring.
    ``evaluable``
        False when the solver returned something that is not a growth rate.
    ``fba_growth_rate``
        The growth rate, or NaN when the row is not evaluable.
    ``predicted``
        "Yes", "No", or None when not evaluable.
    ``category``
        One of true_positive, true_negative, false_positive, false_negative,
        unsure, invalid_solve, excluded. This is the single field downstream
        code should switch on.
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
    ``no_uptake_route``
        True when this row's no-growth prediction rests on a compound having no
        exchange reaction. A sub-count of the negative predictions, not a
        category: these rows are scored like any other.
    ``large_flux_reactions``
        Only present when ``flux_limit`` is given. A frozenset of reaction ids
        that carried more than ``flux_limit`` in this condition's parsimonious
        solution; empty for conditions that did not grow. Take the union across
        rows to get the model-wide count -- see
        :func:`count_large_flux_reactions`.
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
            # Inside the same context manager so the pFBA solve sees exactly
            # the medium the verdict was based on.
            if flux_limit is not None and valid and rate > growth_threshold:
                large_flux = _large_flux_reactions(model, flux_limit)
            else:
                large_flux = frozenset()

        if not valid:
            raw_category = "invalid_solve"
            predicted = None
            predicted_growth = None
        else:
            predicted_growth = rate > growth_threshold
            predicted = "Yes" if predicted_growth else "No"
            raw_category = _classify(row["growth"], predicted_growth)

        excluded = bool(row.get(EXCLUSION_COLUMN, ""))
        record = {
            "missing_exchanges": ", ".join(missing),
            "evaluable": valid,
            "fba_growth_rate": rate if valid else float("nan"),
            "predicted": predicted,
            "raw_category": raw_category,
            "category": "excluded" if excluded else raw_category,
            "discordant": (not excluded) and raw_category in DISCORDANT,
            "near_threshold": bool(
                valid and 0 < rate <= NEAR_THRESHOLD_FACTOR * growth_threshold
            ),
            # Reportable, not a category: this row IS scored, but the no-growth
            # call rests on the model having no way to take the compound up.
            # Note it is False for a row that grows anyway on its other
            # metabolites -- there the missing exchange did not decide anything.
            "no_uptake_route": bool(missing) and predicted_growth is False,
        }
        if flux_limit is not None:
            record["large_flux_reactions"] = large_flux
        records.append(record)

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
        # Not part of the partition: a strict subset of the scored rows, all of
        # them inside true_negative or false_negative. Reported so a caption can
        # say how many negative predictions come from an absent transporter
        # rather than from network structure. Restricted to scored rows on
        # purpose -- an unsure or excluded condition with no uptake route is not
        # a prediction, so including it would make this not add up against the
        # confusion matrix.
        "n_no_uptake_route": int(
            (
                results.get("no_uptake_route", pd.Series(False, index=results.index))
                & results["category"].isin(CONCORDANT + DISCORDANT)
            ).sum()
        ),
        "n_invalid_solve": counts.get("invalid_solve", 0),
        "n_unsure": counts.get("unsure", 0),
        "n_excluded": counts.get("excluded", 0),
        "n_scored": true_positive + true_negative + false_positive + false_negative,
        "true_positive": true_positive,
        "true_negative": true_negative,
        "false_positive": false_positive,
        "false_negative": false_negative,
        # The number plotted as "growth phenotypes matching experimental data".
        # Defined here rather than recomputed by each caller so that the figure,
        # the CI report and the time series cannot disagree about what a match
        # is. It counts only scored rows: an unsure observation, an excluded
        # one, or an unusable solve is not a match either way.
        "matches": true_positive + true_negative,
        "sensitivity": true_positive / observed_growth if observed_growth else np.nan,
        "specificity": (
            true_negative / observed_no_growth if observed_no_growth else np.nan
        ),
    }


def count_interpretable(phenotypes: pd.DataFrame | None = None) -> int:
    """How many conditions could, in principle, be scored against any model.

    A property of the *data*, not of a model: rows with a definite Yes/No
    observation and no exclusion reason. Deliberately not derived from an
    evaluation result, because ``n_scored`` there depends on which exchange
    reactions the model happens to have, and a denominator that grows as the
    model grows cannot be used to compare two models.

    This is the denominator to use for "fraction of growth phenotypes
    reproduced" in the manuscript, and the one plotted in figure 2B.
    """
    if phenotypes is None:
        phenotypes = load_phenotypes()
    definite = phenotypes["growth"].isin(("Yes", "No"))
    not_excluded = phenotypes[EXCLUSION_COLUMN] == ""
    return int((definite & not_excluded).sum())


def count_large_flux_reactions(results: pd.DataFrame) -> int:
    """How many distinct reactions carried a large flux in *any* condition.

    Counts the union rather than summing per condition: one loop that fires in
    every medium is one problem to fix, not sixty.

    Requires ``evaluate_phenotypes(..., flux_limit=...)``.
    """
    if "large_flux_reactions" not in results.columns:
        raise KeyError(
            "results has no 'large_flux_reactions' column; call "
            "evaluate_phenotypes(model, flux_limit=...) to collect it"
        )
    union: set = set()
    for reaction_ids in results["large_flux_reactions"]:
        union.update(reaction_ids)
    return len(union)


def format_summary(summary: dict) -> str:
    """One-paragraph text version of :func:`summarise`, for logs and captions."""
    return (
        f"{summary['n_conditions']} conditions; "
        f"{summary['n_invalid_solve']} invalid solves, "
        f"{summary['n_unsure']} experimentally unsure, "
        f"{summary['n_excluded']} excluded; "
        f"{summary['n_scored']} scored "
        f"({summary['n_no_uptake_route']} of them predicted no-growth because "
        f"the model has no uptake route). "
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

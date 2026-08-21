"""Growth phenotype tests.

Two things are checked here. The model must not grow without a carbon source,
and its predictions across ``data/known_growth_phenotypes.tsv`` must not have
got worse.

The second one used to fail on any mismatch at all. That cannot work for a
model under active curation, where some conditions are always wrong and the
interesting question is whether the set of wrong ones changed. So the currently
accepted mismatches live in ``test_files/expected_phenotype_mismatches.tsv``,
and the test compares against that baseline. A new mismatch fails. A baseline
entry that has started passing also fails, which is the half people leave out —
without it an accidental fix goes unnoticed and the baseline slowly becomes a
list of things that are no longer true.

The simulation itself is not implemented here. It lives in ``tools.phenotypes``
and is shared with ``scripts/generate_growth_report.py``, so the test and the
report cannot disagree about the same model. They previously did: this file set
every uptake to a flat 1000 while the report divided a fixed carbon budget by
the carbon count, so the same condition could pass one and fail the other.
"""

import os
import unittest

import cobra
from gem_utilities import media

from tools.media import MEDIA
from tools.phenotypes import (
    CATEGORIES,
    EXCLUSION_COLUMN,
    compare_to_baseline,
    count_interpretable,
    evaluate_phenotypes,
    format_summary,
    load_expected_mismatches,
    load_phenotypes,
    summarise,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(REPO_ROOT, "model.xml")

#: Loaded once. Reading and parsing the SBML is far slower than solving it.
_MODEL = None


def _model():
    global _MODEL
    if _MODEL is None:
        _MODEL = cobra.io.read_sbml_model(MODEL_PATH)
    return _MODEL


class TestGrowthWithoutCarbon(unittest.TestCase):
    """The model must not grow on a medium with no carbon source.

    A leak here invalidates every other phenotype result, because a model that
    can grow on nothing will appear to grow on anything.
    """

    def test_no_growth_without_carbon(self):
        model = _model()
        with model:
            model.medium = media.clean_media(model, MEDIA["minimal"])
            solution = model.optimize()
            self.assertLessEqual(
                solution.objective_value or 0.0,
                1e-6,
                "the model grows with no carbon source available",
            )


class TestExpectedGrowthPhenotypes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = evaluate_phenotypes(_model())
        cls.summary = summarise(cls.results)
        cls.diff = compare_to_baseline(cls.results, load_expected_mismatches())

    def test_no_new_mismatches(self):
        new = self.diff["new"]
        if not new:
            return
        lines = [
            "",
            f"{len(new)} condition(s) disagree with experiment and are not in "
            f"the accepted baseline:",
            "",
        ]
        lines += [f"  {condition}  ->  {category}" for condition, category in new]
        lines += [
            "",
            "Either fix the model, or -- if this is an accepted consequence of "
            "a deliberate change -- record it with:",
            "    python scripts/update_phenotype_baseline.py",
            "",
            format_summary(self.summary),
        ]
        self.fail("\n".join(lines))

    def test_no_stale_baseline_entries(self):
        """Mismatches that now agree must be removed from the baseline.

        Failing here is good news: the model improved. The test fails anyway so
        that the baseline gets updated, rather than quietly carrying permission
        for a failure that no longer happens.
        """
        resolved = self.diff["resolved"]
        if not resolved:
            return
        lines = ["", f"{len(resolved)} baseline entr(y/ies) no longer mismatch:", ""]
        for condition, category, still_exists in resolved:
            why = "now agrees" if still_exists else "condition no longer in the TSV"
            lines.append(f"  {condition}  (was {category}) -- {why}")
        lines += [
            "",
            "Refresh the baseline to record the improvement:",
            "    python scripts/update_phenotype_baseline.py",
        ]
        self.fail("\n".join(lines))

    def test_mismatch_categories_are_unchanged(self):
        """A false positive turning into a false negative is a real change.

        The condition still disagrees either way, so a test keyed only on
        which conditions mismatch would miss it -- but the two need opposite
        fixes, so it should not pass silently.
        """
        changed = self.diff["changed"]
        if not changed:
            return
        lines = ["", "Mismatch category changed:", ""]
        lines += [
            f"  {condition}: {was} -> {now}" for condition, was, now in changed
        ]
        lines += ["", "    python scripts/update_phenotype_baseline.py"]
        self.fail("\n".join(lines))

    def test_excluded_rows_are_not_scored(self):
        """Rows carrying an exclusion reason must not reach the scoring."""
        scored = self.results[
            self.results["category"].isin(
                ["true_positive", "true_negative", "false_positive", "false_negative"]
            )
        ]
        leaked = sorted(scored[scored[EXCLUSION_COLUMN] != ""]["condition"])
        self.assertFalse(leaked, f"excluded conditions were scored: {leaked}")

    def test_every_phenotype_was_evaluated(self):
        """Guard against a row being silently dropped between file and result."""
        self.assertEqual(
            len(self.results),
            len(load_phenotypes()),
            "evaluate_phenotypes returned a different number of rows than the "
            "phenotype table contains",
        )


class TestSummaryArithmetic(unittest.TestCase):
    """The reported counts must add up.

    These exist because the confusion matrices in figure 2 of the manuscript
    could not be reconciled with the phenotype file by hand. Every number a
    caption quotes comes from :func:`summarise`, so the invariants that make
    those numbers mean what they say are checked here rather than trusted.
    """

    @classmethod
    def setUpClass(cls):
        cls.results = evaluate_phenotypes(_model())
        cls.summary = summarise(cls.results)

    def test_every_row_gets_a_known_category(self):
        unknown = sorted(set(self.results["category"]) - set(CATEGORIES))
        self.assertFalse(
            unknown,
            f"category values not in CATEGORIES: {unknown}. Anything switching "
            f"on category would silently ignore these rows.",
        )

    def test_categories_partition_the_table(self):
        """No row counted twice, none left out."""
        counted = (
            self.summary["n_scored"]
            + self.summary["n_unsure"]
            + self.summary["n_excluded"]
            + self.summary["n_invalid_solve"]
        )
        self.assertEqual(
            counted,
            self.summary["n_conditions"],
            "the category counts do not sum to the number of conditions, so at "
            "least one row is double counted or dropped",
        )

    def test_no_uptake_route_is_a_subset_of_the_negative_predictions(self):
        """The reported count must be reconcilable with the confusion matrix.

        These rows are scored, so each one sits inside a true negative or a
        false negative. If the count exceeds TN + FN it is picking up unsure or
        excluded rows and cannot be quoted alongside the matrix.
        """
        self.assertLessEqual(
            self.summary["n_no_uptake_route"],
            self.summary["true_negative"] + self.summary["false_negative"],
        )

    def test_a_missing_exchange_does_not_by_itself_decide_the_verdict(self):
        """A condition can grow on the compounds it does have.

        ``marine_broth_wo_yeast_and_peptone | Methionine, Pyruvate`` has no
        methionine exchange but grows on the pyruvate, and it grows
        experimentally too. Scoring it off the missing exchange would record a
        false negative for a condition the model gets right, so any row with a
        missing exchange that still grows must not be flagged as a no-uptake
        prediction.
        """
        grew_anyway = self.results[
            (self.results["missing_exchanges"] != "")
            & (self.results["predicted"] == "Yes")
        ]
        self.assertFalse(
            grew_anyway["no_uptake_route"].any(),
            "a condition that grew was flagged as predicting no-growth from a "
            "missing uptake route",
        )

    def test_confusion_matrix_sums_to_the_scored_rows(self):
        self.assertEqual(
            self.summary["true_positive"]
            + self.summary["true_negative"]
            + self.summary["false_positive"]
            + self.summary["false_negative"],
            self.summary["n_scored"],
        )

    def test_matches_are_the_concordant_cells(self):
        """"Matches" is TP + TN and nothing else.

        In particular an unsure observation or an excluded condition is not a
        match. The old time-series code counted excluded conditions, which is
        one of the reasons the figure's numbers could not be reproduced.
        """
        self.assertEqual(
            self.summary["matches"],
            self.summary["true_positive"] + self.summary["true_negative"],
        )

    def test_matches_never_exceed_the_interpretable_denominator(self):
        self.assertLessEqual(
            self.summary["matches"],
            count_interpretable(),
            "more matches than there are interpretable conditions, so the "
            "denominator plotted in figure 2B is wrong",
        )

    def test_excluded_rows_keep_their_underlying_verdict(self):
        """Excluding a row must not destroy the comparison it represents.

        ``raw_category`` is what the row would have scored; only ``category``
        is overwritten. Losing the raw verdict would make it impossible to say
        later whether excluding the row changed the picture.
        """
        excluded = self.results[self.results["category"] == "excluded"]
        self.assertFalse(
            (excluded["raw_category"] == "excluded").any(),
            "an excluded row lost its underlying verdict",
        )


if __name__ == "__main__":
    unittest.main()

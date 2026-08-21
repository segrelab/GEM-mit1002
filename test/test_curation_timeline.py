"""Consistency checks on the committed curation time series.

``curation_process/phenotype_confusion_over_time.csv`` is the file figure 2 of
the manuscript is read off: panel B plots its match counts and panel A quotes
two of its rows as confusion matrices. Nothing recomputes it in CI -- producing
it needs the GitHub API and a solve per PR -- so what is checked here is that
the committed file is internally coherent and was produced by one scoring
definition.

That is the failure this catches. The numbers in the figure originally could not
be reconciled with the phenotype table by hand, and the reason was a file
holding rows scored under two different definitions of "match", with no record
of which was which.

No model is loaded, so this runs in milliseconds.
"""

import os
import unittest

import pandas as pd

from tools.phenotypes import count_interpretable, load_phenotypes

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFUSION_CSV = os.path.join(
    REPO_ROOT, "curation_process", "phenotype_confusion_over_time.csv"
)
SUMMARY_CSV = os.path.join(REPO_ROOT, "curation_process", "growth_match_summary.csv")

#: Kept in sync with ``run_tests_on_prs.SCORING_VERSION`` by
#: :meth:`TestConfusionTimeline.test_scoring_version_matches_the_script`, rather
#: than imported, so that this file does not need cobra installed.
EXPECTED_SCORING_VERSION = 3


def _load(path):
    table = pd.read_csv(path)
    # Rows whose evaluation failed are recorded as ERROR on purpose; they carry
    # no counts to check, and dropping them here is not hiding anything because
    # test_no_errored_rows reports them separately.
    return table[table["Matches"].astype(str) != "ERROR"].copy()


@unittest.skipUnless(
    os.path.exists(CONFUSION_CSV),
    "phenotype_confusion_over_time.csv has not been generated yet; run "
    "curation_process/run_tests_on_prs.py",
)
class TestConfusionTimeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = pd.read_csv(CONFUSION_CSV)
        cls.data = _load(CONFUSION_CSV)
        cls.phenotypes = load_phenotypes()

    def test_one_scoring_version_throughout(self):
        versions = sorted(self.raw["Scoring Version"].dropna().unique())
        self.assertEqual(
            versions,
            [EXPECTED_SCORING_VERSION],
            "the time series mixes scoring definitions (or predates the "
            "version stamp), so its points are not comparable with each other. "
            "Re-run curation_process/run_tests_on_prs.py.",
        )

    def test_scoring_version_matches_the_script(self):
        """Catch the stamp and the constant drifting apart."""
        script = os.path.join(
            REPO_ROOT, "curation_process", "run_tests_on_prs.py"
        )
        with open(script, encoding="utf-8") as handle:
            source = handle.read()
        self.assertIn(
            f"SCORING_VERSION = {EXPECTED_SCORING_VERSION}",
            source,
            "run_tests_on_prs.py declares a different SCORING_VERSION than this "
            "test expects; bump EXPECTED_SCORING_VERSION and regenerate the CSV",
        )

    def test_confusion_matrix_sums_to_scored(self):
        cells = (
            self.data["True Positive"]
            + self.data["True Negative"]
            + self.data["False Positive"]
            + self.data["False Negative"]
        )
        bad = self.data.loc[cells != self.data["Scored"], "PR Number"].tolist()
        self.assertFalse(bad, f"TP+TN+FP+FN != Scored for PR(s) {bad}")

    def test_matches_equal_the_concordant_cells(self):
        concordant = self.data["True Positive"] + self.data["True Negative"]
        bad = self.data.loc[
            concordant != self.data["Matches"], "PR Number"
        ].tolist()
        self.assertFalse(bad, f"Matches != TP+TN for PR(s) {bad}")

    def test_all_categories_sum_to_the_condition_count(self):
        """Every condition lands in exactly one bucket, for every PR."""
        total = (
            self.data["Scored"]
            + self.data["Unsure"]
            + self.data["Excluded"]
            + self.data["Invalid Solve"]
        )
        bad = self.data.loc[total != self.data["Conditions"], "PR Number"].tolist()
        self.assertFalse(
            bad,
            f"the category counts do not add up to Conditions for PR(s) {bad}; "
            f"a condition is being double counted or dropped",
        )

    def test_no_uptake_route_is_within_the_negative_predictions(self):
        """It is a subset of Scored, not another bucket.

        Quoting it next to the confusion matrix only works if it fits inside
        TN + FN. A count larger than that means unsure or excluded rows leaked
        into it.
        """
        negatives = self.data["True Negative"] + self.data["False Negative"]
        bad = self.data.loc[
            self.data["No Uptake Route"] > negatives, "PR Number"
        ].tolist()
        self.assertFalse(
            bad, f"No Uptake Route exceeds TN+FN for PR(s) {bad}"
        )

    def test_condition_count_matches_the_phenotype_table(self):
        """A stale series scored against a different version of the TSV.

        Rows are added to ``known_growth_phenotypes.tsv`` as experiments come
        in. Because the series is cached per PR, an old row can have been
        scored against a smaller table -- which makes the line's early points
        incomparable with its late ones for a reason that has nothing to do
        with the model.
        """
        expected = len(self.phenotypes)
        bad = self.data.loc[
            self.data["Conditions"] != expected, "PR Number"
        ].tolist()
        self.assertFalse(
            bad,
            f"PR(s) {bad} were scored against a different number of conditions "
            f"than the {expected} now in data/known_growth_phenotypes.tsv. "
            f"Re-run run_tests_on_prs.py with FORCE_RERUN = True.",
        )

    def test_denominator_is_constant_and_correct(self):
        expected = count_interpretable(self.phenotypes)
        found = sorted(self.data["Interpretable"].unique())
        self.assertEqual(
            found,
            [expected],
            f"the match denominator should be {expected} for every PR (the "
            f"conditions with a definite Yes/No and no exclusion reason); "
            f"found {found}",
        )

    def test_matches_within_range(self):
        bad = self.data.loc[
            (self.data["Matches"] < 0)
            | (self.data["Matches"] > self.data["Interpretable"]),
            "PR Number",
        ].tolist()
        self.assertFalse(bad, f"Matches outside 0..Interpretable for PR(s) {bad}")

    def test_pr_numbers_are_unique(self):
        counts = self.raw["PR Number"].value_counts()
        self.assertFalse(
            sorted(counts[counts > 1].index),
            f"duplicate PR rows: {sorted(counts[counts > 1].index)}",
        )

    def test_no_errored_rows(self):
        """An ERROR row is a gap in the figure, not a result.

        Not fatal to the analysis, but it should be visible rather than sitting
        in the file unnoticed -- re-running the script retries them.
        """
        errored = self.raw.loc[
            self.raw["Matches"].astype(str) == "ERROR", "PR Number"
        ].tolist()
        self.assertFalse(
            errored,
            f"PR(s) {errored} failed to evaluate and are missing from the "
            f"series; re-run curation_process/run_tests_on_prs.py",
        )

    def test_summary_view_agrees_with_the_full_record(self):
        """``growth_match_summary.csv`` is a view, so it must not disagree."""
        if not os.path.exists(SUMMARY_CSV):
            self.skipTest("growth_match_summary.csv not generated yet")
        summary = _load(SUMMARY_CSV).set_index("PR Number")
        full = self.data.set_index("PR Number")
        shared = summary.index.intersection(full.index)
        self.assertEqual(
            len(shared),
            len(full),
            "the two files cover different sets of PRs; regenerate both by "
            "running run_tests_on_prs.py",
        )
        for column, source in (("Matches", "Matches"), ("Total", "Interpretable")):
            mismatched = shared[
                summary.loc[shared, column].values != full.loc[shared, source].values
            ].tolist()
            self.assertFalse(
                mismatched,
                f"{column} in growth_match_summary.csv disagrees with "
                f"{source} in phenotype_confusion_over_time.csv for PR(s) "
                f"{mismatched}",
            )


if __name__ == "__main__":
    unittest.main()

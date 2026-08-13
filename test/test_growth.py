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
    EXCLUSION_COLUMN,
    compare_to_baseline,
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


if __name__ == "__main__":
    unittest.main()

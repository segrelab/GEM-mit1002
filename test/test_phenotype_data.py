"""Schema tests for data/known_growth_phenotypes.tsv and its baseline.

Deliberately separate from ``test_growth.py``: none of this needs the model, so
it runs in milliseconds and still tells you which file is wrong when the slow
test fails. These guard the parts that other code assumes -- the exclusion
vocabulary, the uniqueness of the condition key, and the baseline pointing at
conditions that actually exist.
"""

import os
import unittest

from tools.phenotypes import (
    CATEGORIES,
    EXCLUSION_COLUMN,
    EXCLUSION_REASONS,
    MISMATCH_COLUMNS,
    condition_key,
    count_interpretable,
    load_expected_mismatches,
    load_phenotypes,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_README = os.path.join(REPO_ROOT, "data", "README.md")


class TestPhenotypeSchema(unittest.TestCase):
    def setUp(self):
        self.phenotypes = load_phenotypes()

    def test_growth_column_vocabulary(self):
        allowed = {"Yes", "No", "Unsure"}
        bad = sorted(set(self.phenotypes["growth"]) - allowed)
        self.assertFalse(
            bad,
            f"growth must be one of {sorted(allowed)}; found {bad}",
        )

    def test_exclude_reason_vocabulary(self):
        used = {r for r in self.phenotypes[EXCLUSION_COLUMN] if r}
        bad = sorted(used - set(EXCLUSION_REASONS))
        self.assertFalse(
            bad,
            f"{EXCLUSION_COLUMN} values not in the controlled vocabulary: {bad}. "
            f"Allowed: {sorted(EXCLUSION_REASONS)}. Adding a category means "
            f"updating tools/phenotypes.py and data/README.md together.",
        )

    def test_every_reason_is_documented(self):
        """The vocabulary and the README have to be changed together."""
        with open(DATA_README, encoding="utf-8") as handle:
            readme = handle.read()
        undocumented = sorted(
            reason for reason in EXCLUSION_REASONS if f"`{reason}`" not in readme
        )
        self.assertFalse(
            undocumented,
            f"these exclusion reasons are not described in data/README.md: "
            f"{undocumented}",
        )

    def test_condition_keys_are_unique(self):
        """``minimal_media`` + ``c_source`` identifies a row.

        The expected-mismatch baseline is keyed on this pair, so a duplicate
        would make a baseline entry ambiguous.
        """
        counts = self.phenotypes["condition"].value_counts()
        duplicated = sorted(counts[counts > 1].index)
        self.assertFalse(
            duplicated,
            f"duplicate (minimal_media, c_source) keys: {duplicated}",
        )

    def test_met_ids_are_present_and_well_formed(self):
        bad = [
            (row["condition"], row["met_id"])
            for _, row in self.phenotypes.iterrows()
            if not row["met_id"]
            or not all(m.startswith("cpd") for m in row["met_id"])
        ]
        self.assertFalse(bad, f"rows with missing or malformed met_id: {bad}")


class TestInterpretableCount(unittest.TestCase):
    """The denominator quoted in the manuscript and plotted in figure 2B.

    Checked here, against the file, because it is a property of the data and
    must not depend on which model it happens to be compared against.
    """

    def setUp(self):
        self.phenotypes = load_phenotypes()

    def test_categories_are_distinct(self):
        self.assertEqual(
            len(CATEGORIES),
            len(set(CATEGORIES)),
            "a category name appears twice, so counts built from CATEGORIES "
            "would double count those rows",
        )

    def test_interpretable_excludes_unsure_and_excluded_rows(self):
        expected = sum(
            1
            for _, row in self.phenotypes.iterrows()
            if row["growth"] in ("Yes", "No") and not row[EXCLUSION_COLUMN]
        )
        self.assertEqual(count_interpretable(self.phenotypes), expected)

    def test_interpretable_is_a_strict_subset_when_rows_are_held_out(self):
        """A sanity check on the direction of the filter.

        If the count ever equalled the row count while the file still contains
        unsure or excluded rows, the filter has stopped filtering.
        """
        held_out = sum(
            1
            for _, row in self.phenotypes.iterrows()
            if row["growth"] not in ("Yes", "No") or row[EXCLUSION_COLUMN]
        )
        if not held_out:
            self.skipTest("no unsure or excluded rows in the table")
        self.assertLess(count_interpretable(self.phenotypes), len(self.phenotypes))


class TestExpectedMismatches(unittest.TestCase):
    def setUp(self):
        self.phenotypes = load_phenotypes()
        self.expected = load_expected_mismatches()

    def test_baseline_has_the_expected_columns(self):
        if self.expected.empty and not len(self.expected.columns):
            self.skipTest("no baseline file yet")
        missing = [c for c in MISMATCH_COLUMNS if c not in self.expected.columns]
        self.assertFalse(missing, f"baseline is missing columns: {missing}")

    def test_baseline_rows_refer_to_real_conditions(self):
        """A baseline entry naming a condition that no longer exists is stale.

        This happens when a phenotype is renamed or removed and the baseline is
        not regenerated, and it would otherwise sit there forever excusing a
        mismatch that cannot occur.
        """
        known = set(self.phenotypes["condition"])
        orphaned = sorted(set(self.expected.get("condition", [])) - known)
        self.assertFalse(
            orphaned,
            f"baseline entries with no matching phenotype row: {orphaned}. "
            f"Regenerate with scripts/update_phenotype_baseline.py.",
        )

    def test_baseline_does_not_list_excluded_conditions(self):
        """Excluded rows are never scored, so they can never be a mismatch."""
        excluded = {
            condition_key(row["minimal_media"], row["c_source"])
            for _, row in self.phenotypes.iterrows()
            if row[EXCLUSION_COLUMN]
        }
        overlap = sorted(set(self.expected.get("condition", [])) & excluded)
        self.assertFalse(
            overlap,
            f"baseline lists conditions that are excluded from scoring, so the "
            f"entries are dead: {overlap}",
        )


if __name__ == "__main__":
    unittest.main()

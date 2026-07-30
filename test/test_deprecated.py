"""Tests for the deprecated-identifier lists in data/deprecatedIdentifiers/.

The test that matters most is :meth:`TestDeprecatedNotInModel` -- it is what
stops a curator re-adding something the project already decided to remove, and
stops the two halves of the record drifting apart. The rest guard the schema so
the files stay machine-readable.
"""

import csv
import os
import re
import unittest

import cobra

from tools.deprecate import (
    COLUMNS,
    METABOLITES_TSV,
    MODEL_PATH,
    NOTES_INFO_KEY,
    NOTES_METABOLITE_KEY,
    NOTES_REACTION_KEY,
    REACTIONS_TSV,
    REASONS,
    REASONS_REQUIRING_REPLACEMENT,
    DeprecationError,
    build_notes_dict,
    read_records,
    stamp_pr_number,
    strip_sbml_prefix,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEPRECATED_README = os.path.join(
    REPO_ROOT, "data", "deprecatedIdentifiers", "README.md"
)


def _load_model():
    return cobra.io.read_sbml_model(MODEL_PATH)


class TestDeprecatedNotInModel(unittest.TestCase):
    """A deprecated identifier must not be back in the model.

    If this fails, either someone re-added a reaction the project deliberately
    removed, or a removal was reverted without updating the TSV. Both are worth
    a conversation, which is the point.
    """

    @classmethod
    def setUpClass(cls):
        cls.model = _load_model()

    def test_deprecated_reactions_absent(self):
        listed = {r.id for r in read_records(REACTIONS_TSV)}
        present = {strip_sbml_prefix(r.id) for r in self.model.reactions}
        resurrected = sorted(listed & present)
        self.assertEqual(
            [],
            resurrected,
            msg=(
                f"{len(resurrected)} reaction(s) are in the model but listed as "
                f"deprecated: {resurrected}. Either remove them again, or delete "
                f"their row from {os.path.relpath(REACTIONS_TSV, REPO_ROOT)} and "
                f"explain in your PR why the earlier decision was wrong."
            ),
        )

    def test_deprecated_metabolites_absent(self):
        listed = {r.id for r in read_records(METABOLITES_TSV)}
        present = {strip_sbml_prefix(m.id) for m in self.model.metabolites}
        resurrected = sorted(listed & present)
        self.assertEqual(
            [],
            resurrected,
            msg=(
                f"{len(resurrected)} metabolite(s) are in the model but listed as "
                f"deprecated: {resurrected}. Either remove them again, or delete "
                f"their row from {os.path.relpath(METABOLITES_TSV, REPO_ROOT)} and "
                f"explain in your PR why the earlier decision was wrong."
            ),
        )


class TestDeprecatedSchema(unittest.TestCase):
    """The TSVs must stay well-formed and machine-readable."""

    def test_headers_exact(self):
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                self.assertTrue(os.path.exists(path), msg=f"{path} is missing")
                with open(path, newline="", encoding="utf-8") as handle:
                    header = next(csv.reader(handle, delimiter="\t"))
                self.assertEqual(
                    COLUMNS,
                    header,
                    msg=(
                        f"{os.path.basename(path)} header drifted from the schema. "
                        f"Expected {COLUMNS}."
                    ),
                )

    def test_no_duplicate_ids(self):
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                ids = [r.id for r in read_records(path)]
                dupes = sorted({i for i in ids if ids.count(i) > 1})
                self.assertEqual(
                    [],
                    dupes,
                    msg=(
                        f"duplicate rows in {os.path.basename(path)}: {dupes}. "
                        f"One row per identifier; if it was removed twice, keep "
                        f"the row and note both PRs."
                    ),
                )

    def test_reasons_in_vocabulary(self):
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                bad = sorted(
                    {
                        (r.id, r.reason)
                        for r in read_records(path)
                        if r.reason not in REASONS
                    }
                )
                self.assertEqual(
                    [],
                    bad,
                    msg=(
                        f"reason values outside the controlled vocabulary in "
                        f"{os.path.basename(path)}: {bad}. Allowed: "
                        f"{', '.join(REASONS)}. To add a category, update REASONS "
                        f"in tools/deprecate.py and the table in the README."
                    ),
                )

    def test_ids_have_no_sbml_prefix(self):
        """Stored bare, so they match what COBRApy reports."""
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                with open(path, newline="", encoding="utf-8") as handle:
                    raw = [row["id"] for row in csv.DictReader(handle, delimiter="\t")]
                prefixed = sorted(i for i in raw if re.match(r"^(R_|M_|G_)", i))
                self.assertEqual(
                    [],
                    prefixed,
                    msg=(
                        f"identifiers in {os.path.basename(path)} carry an SBML "
                        f"prefix: {prefixed}. Store them bare (rxn00196_c0, not "
                        f"R_rxn00196_c0)."
                    ),
                )

    def test_dates_are_iso(self):
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                bad = sorted(
                    (r.id, r.date)
                    for r in read_records(path)
                    if r.date and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", r.date)
                )
                self.assertEqual(
                    [], bad, msg=f"dates must be YYYY-MM-DD: {bad}"
                )

    def test_pr_format(self):
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                bad = sorted(
                    (r.id, r.pr)
                    for r in read_records(path)
                    if r.pr and not re.fullmatch(r"#\d+", r.pr)
                )
                self.assertEqual(
                    [], bad, msg=f"pr must look like '#123': {bad}"
                )

    def test_notes_are_single_line(self):
        """Long-form reasoning belongs in the PR, not in the TSV."""
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                bad = sorted(
                    r.id for r in read_records(path) if len(r.notes) > 200
                )
                self.assertEqual(
                    [],
                    bad,
                    msg=(
                        f"notes longer than 200 characters in "
                        f"{os.path.basename(path)}: {bad}. Summarise here and put "
                        f"the detail in the PR."
                    ),
                )


class TestReplacedBy(unittest.TestCase):
    """``replaced_by`` is the most actionable column, so it gets checked."""

    @classmethod
    def setUpClass(cls):
        cls.model = _load_model()
        cls.reaction_ids = {strip_sbml_prefix(r.id) for r in cls.model.reactions}
        cls.metabolite_ids = {strip_sbml_prefix(m.id) for m in cls.model.metabolites}

    @staticmethod
    def _targets(record):
        return [t.strip() for t in record.replaced_by.split(";") if t.strip()]

    def test_required_when_reason_is_relative(self):
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                bad = sorted(
                    (r.id, r.reason)
                    for r in read_records(path)
                    if r.reason in REASONS_REQUIRING_REPLACEMENT and not r.replaced_by
                )
                self.assertEqual(
                    [],
                    bad,
                    msg=(
                        f"reason in {REASONS_REQUIRING_REPLACEMENT} is a claim about "
                        f"another identifier, so replaced_by must be filled in: {bad}"
                    ),
                )

    def test_not_self_referential(self):
        for path in (REACTIONS_TSV, METABOLITES_TSV):
            with self.subTest(path=os.path.basename(path)):
                bad = sorted(
                    r.id for r in read_records(path) if r.id in self._targets(r)
                )
                self.assertEqual(
                    [], bad, msg=f"replaced_by points at itself: {bad}"
                )

    def test_targets_resolve(self):
        """A replacement should be in the model, or itself deprecated.

        A pointer to something that never existed is worse than no pointer, so
        this catches typos. Chains (A replaced by B, B later replaced by C) are
        allowed -- the reader can follow them.
        """
        for path, in_model in (
            (REACTIONS_TSV, self.reaction_ids),
            (METABOLITES_TSV, self.metabolite_ids),
        ):
            with self.subTest(path=os.path.basename(path)):
                records = read_records(path)
                deprecated = {r.id for r in records}
                dangling = sorted(
                    {
                        (r.id, t)
                        for r in records
                        for t in self._targets(r)
                        if t not in in_model and t not in deprecated
                    }
                )
                self.assertEqual(
                    [],
                    dangling,
                    msg=(
                        f"replaced_by targets that are neither in the model nor "
                        f"themselves deprecated (likely typos): {dangling}"
                    ),
                )


class TestNotesMirror(unittest.TestCase):
    """The model's <notes> mirror must agree with the TSVs.

    ``scripts/export_model.py`` regenerates it, so a failure here means the
    model file was edited without re-running the export.
    """

    @classmethod
    def setUpClass(cls):
        cls.model = _load_model()
        cls.expected = build_notes_dict()

    def _listed(self, key):
        raw = self.model.notes.get(key, "")
        return {part.strip() for part in raw.split(";") if part.strip()}

    def test_mirror_matches_tsv(self):
        if not self.expected:
            self.skipTest("no deprecated identifiers recorded yet")
        for key in (NOTES_REACTION_KEY, NOTES_METABOLITE_KEY):
            if key not in self.expected:
                continue
            with self.subTest(key=key):
                expected = {
                    p.strip() for p in self.expected[key].split(";") if p.strip()
                }
                actual = self._listed(key)
                self.assertEqual(
                    expected,
                    actual,
                    msg=(
                        f"model notes {key} disagrees with the TSV. Missing from "
                        f"the model: {sorted(expected - actual)}; unexpectedly "
                        f"present: {sorted(actual - expected)}. Run "
                        f"`python -m tools.deprecate sync-notes` to regenerate."
                    ),
                )

    def test_pointer_present(self):
        if not self.expected:
            self.skipTest("no deprecated identifiers recorded yet")
        self.assertIn(
            NOTES_INFO_KEY,
            self.model.notes,
            msg=(
                "the model notes should point at the full table, so someone with "
                "only model.xml can find the reasons. Run "
                "`python -m tools.deprecate sync-notes`."
            ),
        )

    def test_mirror_survives_cobrapy_round_trip(self):
        """The reason we use <notes> rather than a custom <annotation>.

        COBRApy drops foreign-namespace annotation blocks on write, so a custom
        section would vanish the first time anyone round-tripped the model --
        including our own export script. Notes survive. This test would catch a
        COBRApy change that broke that assumption.
        """
        if not self.expected:
            self.skipTest("no deprecated identifiers recorded yet")
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "round_trip.xml")
            cobra.io.write_sbml_model(self.model, out)
            reloaded = cobra.io.read_sbml_model(out)
        for key in (NOTES_REACTION_KEY, NOTES_METABOLITE_KEY, NOTES_INFO_KEY):
            if key not in self.expected:
                continue
            with self.subTest(key=key):
                self.assertEqual(
                    self.model.notes.get(key),
                    reloaded.notes.get(key),
                    msg=(
                        f"{key} did not survive a COBRApy write/read cycle. If "
                        f"COBRApy changed how it handles model notes, the mirroring "
                        f"strategy in data/deprecatedIdentifiers/README.md needs "
                        f"revisiting."
                    ),
                )


class TestStampPrNumber(unittest.TestCase):
    """CI fills in the pr column, because you cannot know it in advance.

    The property that matters is that stamping only fills blank cells: a row
    deliberately pointing at a different PR, or at an issue, must survive.
    """

    HEADER = "\t".join(COLUMNS) + "\n"

    def _tsv(self, tmp, name, rows):
        path = os.path.join(tmp, name)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(self.HEADER)
            for row in rows:
                handle.write("\t".join(row) + "\n")
        return path

    def test_fills_blanks_and_preserves_existing(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            rxns = self._tsv(
                tmp,
                "deprecatedReactions.tsv",
                [
                    ("rxnAAA_c0", "blank", "dead_end", "", "", "2026-01-01", ""),
                    ("rxnBBB_c0", "issue", "dead_end", "", "#316", "2026-01-01", ""),
                ],
            )
            mets = self._tsv(tmp, "deprecatedMetabolites.tsv", [])
            stamped = stamp_pr_number("317", reactions_tsv=rxns, metabolites_tsv=mets)

            self.assertEqual({rxns: ["rxnAAA_c0"]}, stamped)
            by_id = {r.id: r for r in read_records(rxns)}
            self.assertEqual("#317", by_id["rxnAAA_c0"].pr)
            self.assertEqual(
                "#316",
                by_id["rxnBBB_c0"].pr,
                msg="stamping must not overwrite a row already attributed elsewhere",
            )

    def test_is_idempotent(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            rxns = self._tsv(
                tmp,
                "deprecatedReactions.tsv",
                [("rxnAAA_c0", "blank", "dead_end", "", "", "2026-01-01", "")],
            )
            mets = self._tsv(tmp, "deprecatedMetabolites.tsv", [])
            stamp_pr_number("317", reactions_tsv=rxns, metabolites_tsv=mets)
            first = open(rxns, encoding="utf-8").read()
            self.assertEqual(
                {}, stamp_pr_number("999", reactions_tsv=rxns, metabolites_tsv=mets)
            )
            self.assertEqual(
                first,
                open(rxns, encoding="utf-8").read(),
                msg="a second CI run must not change anything",
            )

    def test_rejects_non_numeric(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            rxns = self._tsv(tmp, "deprecatedReactions.tsv", [])
            mets = self._tsv(tmp, "deprecatedMetabolites.tsv", [])
            for bad in ("", "not-a-number", "#abc"):
                with self.subTest(value=bad):
                    with self.assertRaises(DeprecationError):
                        stamp_pr_number(
                            bad, reactions_tsv=rxns, metabolites_tsv=mets
                        )


class TestVocabularyDocumented(unittest.TestCase):
    """Every allowed reason must be explained to contributors."""

    def test_every_reason_in_readme(self):
        self.assertTrue(
            os.path.exists(DEPRECATED_README), msg=f"{DEPRECATED_README} is missing"
        )
        text = open(DEPRECATED_README, encoding="utf-8").read()
        undocumented = sorted(r for r in REASONS if f"`{r}`" not in text)
        self.assertEqual(
            [],
            undocumented,
            msg=(
                f"these reason values are allowed by tools/deprecate.py but not "
                f"documented in data/deprecatedIdentifiers/README.md: "
                f"{undocumented}. A vocabulary nobody can look up is not a "
                f"vocabulary."
            ),
        )


if __name__ == "__main__":
    unittest.main()

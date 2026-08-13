"""Record the current phenotype mismatches as the accepted baseline.

``test/test_growth.py`` compares the model against
``test/test_files/expected_phenotype_mismatches.tsv`` rather than demanding a
perfect match, because a model under curation always has some conditions wrong
and the useful question is whether that set changed.

Run this when the set has changed for a reason you accept -- you fixed
something, or you made a deliberate trade -- then commit the result alongside
the change that caused it. Reviewing that diff is the point: it is the record
of what a curation step did to the model's agreement with experiment.

Never run this from CI. Automatically accepting whatever the model currently
does would turn the test into a rubber stamp.

    python scripts/update_phenotype_baseline.py
    python scripts/update_phenotype_baseline.py --dry-run
"""

import argparse
import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from tools.phenotypes import (  # noqa: E402
    EXPECTED_MISMATCHES_TSV,
    compare_to_baseline,
    evaluate_phenotypes,
    format_summary,
    load_expected_mismatches,
    summarise,
    write_expected_mismatches,
)

MODEL_PATH = os.path.join(PROJECT_ROOT, "model.xml")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="python scripts/update_phenotype_baseline.py",
        description=__doc__.split("\n\n")[0],
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show what would change and write nothing",
    )
    parser.add_argument(
        "--notes",
        default="",
        help="one short line recorded against every entry, e.g. the PR subject",
    )
    parser.add_argument("--model", default=MODEL_PATH)
    args = parser.parse_args(argv)

    import cobra

    model = cobra.io.read_sbml_model(args.model)
    results = evaluate_phenotypes(model)
    summary = summarise(results)
    diff = compare_to_baseline(results, load_expected_mismatches())

    print(format_summary(summary))
    print()

    for label, rows in (
        ("newly mismatching", [(c, cat) for c, cat in diff["new"]]),
        ("no longer mismatching", [(c, cat) for c, cat, _ in diff["resolved"]]),
        (
            "changed category",
            [(c, f"{was} -> {now}") for c, was, now in diff["changed"]],
        ),
    ):
        print(f"{label}: {len(rows)}")
        for condition, detail in rows:
            print(f"    {condition}  {detail}")

    if not any(diff.values()):
        print("\nBaseline already matches the current model; nothing to write.")
        return 0

    if args.dry_run:
        print("\n--dry-run: nothing written.")
        return 0

    written = write_expected_mismatches(
        results, EXPECTED_MISMATCHES_TSV, notes=args.notes
    )
    rel = os.path.relpath(EXPECTED_MISMATCHES_TSV, PROJECT_ROOT)
    print(f"\nWrote {len(written)} accepted mismatch(es) to {rel}")
    print("Commit this alongside the change that caused it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

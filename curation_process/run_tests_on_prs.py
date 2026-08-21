"""Replay the growth-phenotype evaluation across every merged PR.

This is what figure 2B of the manuscript plots: how the model's agreement with
the known growth phenotypes moved as curation proceeded. For each merged PR
that touched ``model.xml``, the model at that PR is fetched from GitHub,
scored, and one row is written per PR.

Two outputs, written from the same in-memory record so they cannot disagree:

``phenotype_confusion_over_time.csv``
    The full record: the confusion matrix (TP/TN/FP/FN) plus every unscored
    category, sensitivity, specificity, and model size. This is the source of
    truth -- the confusion matrices in figure 2A should be read off this file
    rather than counted by hand.

``growth_match_summary.csv``
    The narrow view the plotting script reads. A subset of the same numbers.

The scoring itself is NOT implemented here. It lives in ``tools.phenotypes``
and is shared with ``test/test_growth.py``, so the time series and the CI test
cannot disagree about what a match is. They previously did, and the old inlined
loop here got three things wrong that mattered for the figure:

1. ``exclude_reason`` was ignored. The six conditions from the failed
   no-nitrogen control (``marine_broth_wo_yeast_and_peptone_no_n``) were scored
   as matches or mismatches even though the experiment cannot be interpreted.
   This is the main reason the plotted totals could not be reconciled with the
   phenotype file by hand.

2. A condition with a missing exchange reaction was scored as a "No" by fiat,
   without looking at what the model did. Those conditions are scored again
   now, but on the model's actual solution -- which is not the same thing,
   because a multi-compound condition can still grow on the compounds that
   *are* present. See the "Missing exchange reactions" section of
   ``tools/phenotypes.py`` for the two rows where this changes the verdict.

3. Results were cached per PR with no record of how they were scored, so a
   change to the scoring rules left a file mixing old and new numbers. Hence
   ``SCORING_VERSION`` below: a stored row scored under an older version is
   re-run automatically.
"""

import json
import os
import subprocess
import sys
import warnings

import cobra
import pandas as pd

# ---- CONFIG: edit these, then run `python run_tests_on_prs.py` ----
# Range of merged PRs to evaluate. PR_END = None goes up to the latest PR.
# For the initial model construction to growing on everything use 89 - 212
# To highlight the acetate/leucine/isoleucine fixes use 273 - 293
# For initial to after fixing the TCA cycle fluxes use 89 - 344
PR_START = 89
PR_END = None
# Reaction ID for the biomass reaction
BIOMASS_RXN_ID = "bio1_biomass"
# False = incremental (only new PRs, PRs that previously errored, and PRs whose
# stored result was scored under an older SCORING_VERSION).
# True = re-run every PR in the range regardless.
FORCE_RERUN = False

#: Bump this whenever a change to the scoring rules makes previously stored
#: rows non-comparable with new ones. Rows stamped with an older version are
#: re-evaluated automatically, so the file is never a mix of two definitions.
#:
#: 1 = the original inlined loop (exclusions ignored; a missing exchange scored
#:     by fiat as a "No", so a multi-compound condition that the model actually
#:     grows on was still recorded as no-growth).
#: 2 = scored by tools.phenotypes, with conditions whose compound has no
#:     exchange reaction held out entirely as ``no_exchange``.
#: 3 = those conditions scored again, but on the model's actual solution rather
#:     than on the missing exchange. An absent transporter that a genome search
#:     found no evidence for is a prediction of no growth, not a gap, and it is
#:     most of what the model has to say about failure to grow. Infeasible
#:     solves (no carbon source the model can take up, against a forced ATP
#:     maintenance bound) are read as no growth rather than as solver failures.
SCORING_VERSION = 3

# FILE PATHS
FILE_PATH = os.path.dirname(os.path.abspath(__file__))
REPO_PATH = os.path.dirname(FILE_PATH)

sys.path.insert(0, str(REPO_PATH))

from tools.phenotypes import (  # noqa: E402
    DEFAULT_FLUX_LIMIT,
    count_interpretable,
    count_large_flux_reactions,
    evaluate_phenotypes,
    load_phenotypes,
    summarise,
)

SUMMARY_FILE = os.path.join(FILE_PATH, "growth_match_summary.csv")
CONFUSION_FILE = os.path.join(FILE_PATH, "phenotype_confusion_over_time.csv")

TEMP_MODEL = os.path.join(REPO_PATH, "temp_model.xml")

#: Column order of the confusion-matrix file. Explicit so the file's schema is
#: reviewable in a diff instead of following whatever order a dict happened to
#: have.
CONFUSION_COLUMNS = [
    "PR Number",
    "Date Opened",
    "Date Merged",
    "Scoring Version",
    "Reactions",
    "Metabolites",
    "Genes",
    # Confusion matrix over the scored conditions.
    "True Positive",
    "True Negative",
    "False Positive",
    "False Negative",
    # Matches = TP + TN. Interpretable is the model-independent denominator,
    # so Matches / Interpretable is comparable between any two rows.
    "Matches",
    "Interpretable",
    "% Match",
    # Where the rest of the conditions went. Scored + Unsure + Excluded +
    # Invalid Solve must equal Conditions; if it does not, something is being
    # double counted.
    "Scored",
    "Conditions",
    "Unsure",
    "Excluded",
    "Invalid Solve",
    # NOT part of that sum: a subset of Scored. How many of the model's
    # no-growth predictions rest on it having no uptake route for the compound,
    # rather than on network structure downstream of uptake.
    "No Uptake Route",
    "Sensitivity",
    "Specificity",
    "Unbounded Flux Reactions",
]

#: Column order of the narrow file the plot reads. "Total" is kept under its
#: old name for compatibility but now holds the interpretable denominator.
SUMMARY_COLUMNS = [
    "PR Number",
    "Date Opened",
    "Date Merged",
    "Scoring Version",
    "Reactions",
    "Metabolites",
    "Genes",
    "Matches",
    "Total",
    "% Match",
    "Unbounded Flux Reactions",
]


def score_current_temp_model(
    phenotypes: pd.DataFrame,
    n_interpretable: int,
    flux_limit: float = DEFAULT_FLUX_LIMIT,
    biomass_rxn_id: str = BIOMASS_RXN_ID,
) -> dict:
    """Score ``temp_model.xml`` and return one flat record of results.

    ``biomass_rxn_id`` is not used to set the objective -- that comes from the
    SBML -- but its presence is checked, so a PR that renamed or dropped the
    biomass reaction fails loudly instead of quietly reporting that the model
    stopped growing on everything.
    """
    model = cobra.io.read_sbml_model(TEMP_MODEL)

    if biomass_rxn_id not in {reaction.id for reaction in model.reactions}:
        raise ValueError(
            f"model has no reaction {biomass_rxn_id!r}; growth results for this "
            f"PR would not be comparable with the rest of the series"
        )

    results = evaluate_phenotypes(model, phenotypes, flux_limit=flux_limit)
    summary = summarise(results)

    return {
        "Scoring Version": SCORING_VERSION,
        "Reactions": len(model.reactions),
        "Metabolites": len(model.metabolites),
        "Genes": len(model.genes),
        "True Positive": summary["true_positive"],
        "True Negative": summary["true_negative"],
        "False Positive": summary["false_positive"],
        "False Negative": summary["false_negative"],
        "Matches": summary["matches"],
        "Interpretable": n_interpretable,
        "% Match": round(100 * summary["matches"] / max(n_interpretable, 1), 2),
        "Scored": summary["n_scored"],
        "Conditions": summary["n_conditions"],
        "No Uptake Route": summary["n_no_uptake_route"],
        "Unsure": summary["n_unsure"],
        "Excluded": summary["n_excluded"],
        "Invalid Solve": summary["n_invalid_solve"],
        "Sensitivity": (
            round(summary["sensitivity"], 4)
            if pd.notna(summary["sensitivity"])
            else ""
        ),
        "Specificity": (
            round(summary["specificity"], 4)
            if pd.notna(summary["specificity"])
            else ""
        ),
        "Unbounded Flux Reactions": count_large_flux_reactions(results),
    }


def _load_existing() -> dict:
    """Stored confusion rows, keyed by PR number. Missing file reads as empty."""
    if not os.path.exists(CONFUSION_FILE):
        return {}
    stored = pd.read_csv(CONFUSION_FILE)
    return {int(row["PR Number"]): row.to_dict() for _, row in stored.iterrows()}


def _needs_rerun(row: dict) -> bool:
    """True when a stored row cannot be reused.

    Either it errored, or it was produced by an older scoring definition and so
    is not comparable with the rows this run would add.
    """
    if str(row.get("Matches")) == "ERROR":
        return True
    try:
        return int(row.get("Scoring Version", 0)) != SCORING_VERSION
    except (TypeError, ValueError):
        return True


def run_tests_on_prs(
    pr_start: int = PR_START,
    pr_end: int | None = PR_END,
    flux_limit: float = DEFAULT_FLUX_LIMIT,
    biomass_rxn_id: str = BIOMASS_RXN_ID,
    force_rerun: bool = False,
    prs_to_skip=(),
):
    """Evaluate merged PRs and write the two summary files.

    Incremental by default: a PR is (re-)evaluated when it has no stored row,
    when its stored row is an ERROR, or when its stored row was scored under an
    older ``SCORING_VERSION``. Results for PRs outside the requested range are
    preserved untouched.
    """
    phenotypes = load_phenotypes()
    n_interpretable = count_interpretable(phenotypes)
    print(
        f"Scoring against {len(phenotypes)} conditions, "
        f"{n_interpretable} of them interpretable "
        f"(definite Yes/No and not excluded)."
    )

    existing_results = _load_existing()
    reusable = {
        pr_number
        for pr_number, row in existing_results.items()
        if not _needs_rerun(row)
    }
    stale = sorted(set(existing_results) - reusable)
    if stale and not force_rerun:
        print(
            f"{len(stale)} stored row(s) errored or were scored under an older "
            f"definition; they will be re-evaluated."
        )

    # Get a list of PRs merged into the dev branch (could change to main, by
    # providing target_branch="main")
    # By only looking at dev we miss a few PRs that were merged into main,
    # (89-117) but this is easier to automate and still captures the majority
    # of changes to the model
    all_pr_entries = get_prs_by_target()
    if pr_end is None:
        pr_end = max(pr_entry["number"] for pr_entry in all_pr_entries)
    pr_entries_to_check = [
        pr_entry
        for pr_entry in all_pr_entries
        if pr_start <= pr_entry["number"] <= pr_end
        and pr_entry["number"] not in prs_to_skip
        and (force_rerun or pr_entry["number"] not in reusable)
    ]
    # Only PRs that changed the model can change the results
    pull_requests = [
        pr_entry
        for pr_entry in pr_entries_to_check
        if is_model_changed_in_pr(pr_entry["number"])
    ]

    if not pull_requests:
        print("No PRs to evaluate; existing results are up to date.")
    else:
        print(f"Evaluating {len(pull_requests)} PR(s)...")

    results_by_pr = dict(existing_results)

    for pr in pull_requests:
        print(f"\n--- Evaluating PR #{pr['number']} ---")
        record = {
            "PR Number": pr["number"],
            "Date Opened": pr["createdAt"],
            "Date Merged": pr["mergedAt"],
        }
        try:
            fetch_model_at_pr(pr["number"])
            record.update(
                score_current_temp_model(
                    phenotypes,
                    n_interpretable,
                    flux_limit=flux_limit,
                    biomass_rxn_id=biomass_rxn_id,
                )
            )
            print(
                f"    matches {record['Matches']}/{record['Interpretable']}"
                f"  (TP {record['True Positive']}, TN {record['True Negative']}, "
                f"FP {record['False Positive']}, FN {record['False Negative']}; "
                f"{record['No Uptake Route']} from an absent uptake route)"
            )
        except Exception as error:
            print(f"Error with PR #{pr['number']}: {error}")
            record.update(
                {
                    column: "ERROR"
                    for column in CONFUSION_COLUMNS
                    if column not in record
                }
            )
            # Keep the version stamp real even on failure, so a retry is driven
            # by the ERROR marker rather than by the row looking stale.
            record["Scoring Version"] = SCORING_VERSION
        results_by_pr[pr["number"]] = record

    if os.path.exists(TEMP_MODEL):
        os.remove(TEMP_MODEL)

    write_results(results_by_pr)


def fetch_model_at_pr(pr_number: int) -> None:
    """Download ``model.xml`` as of ``pr_number`` into ``temp_model.xml``.

    ``subprocess.run`` is checked here. It was not, so a failed download left
    the *previous* PR's model on disk and the run silently scored the same
    model twice.
    """
    with open(TEMP_MODEL, "w") as handle:
        subprocess.run(
            [
                "gh",
                "api",
                f"repos/:owner/:repo/contents/model.xml?ref=pull/{pr_number}/head",
                "-H",
                "Accept: application/vnd.github.v3.raw",
            ],
            stdout=handle,
            check=True,
        )
    if os.path.getsize(TEMP_MODEL) == 0:
        raise RuntimeError(f"downloaded model.xml for PR #{pr_number} is empty")


#: Columns that hold counts. Written as integers rather than letting pandas
#: promote them to float, which it does as soon as one row is missing a value:
#: "1426.0 reactions" in a file that feeds figure captions is a nuisance.
INTEGER_COLUMNS = [
    "PR Number",
    "Scoring Version",
    "Reactions",
    "Metabolites",
    "Genes",
    "True Positive",
    "True Negative",
    "False Positive",
    "False Negative",
    "Matches",
    "Interpretable",
    "Scored",
    "Conditions",
    "No Exchange",
    "Unsure",
    "Excluded",
    "Invalid Solve",
    "No Uptake Route",
    "Unbounded Flux Reactions",
    "Total",
]


def _as_nullable_ints(table: pd.DataFrame) -> pd.DataFrame:
    """Cast count columns to nullable ints where every value permits it.

    A column holding "ERROR" for a failed PR is left alone rather than
    coerced, because turning a recorded failure into a blank would hide it.
    """
    for column in INTEGER_COLUMNS:
        if column not in table.columns:
            continue
        try:
            table[column] = table[column].astype("Int64")
        except (TypeError, ValueError):
            pass
    return table


def write_results(results_by_pr: dict) -> None:
    """Write both output files from one set of records, sorted by PR number."""
    records = [results_by_pr[number] for number in sorted(results_by_pr)]
    confusion = _as_nullable_ints(pd.DataFrame(records).reindex(columns=CONFUSION_COLUMNS))
    confusion.to_csv(CONFUSION_FILE, index=False, lineterminator="\n")

    # The narrow file is a view of the same numbers, not a second measurement.
    summary = confusion.copy()
    summary["Total"] = summary["Interpretable"]
    summary = _as_nullable_ints(summary.reindex(columns=SUMMARY_COLUMNS))
    summary.to_csv(SUMMARY_FILE, index=False, lineterminator="\n")

    print(f"\nWrote {len(records)} row(s) to:")
    print(f"  {CONFUSION_FILE}")
    print(f"  {SUMMARY_FILE}")


def get_prs_by_target(target_branch="dev"):
    # This requires the GitHub CLI 'gh' to be installed and authenticated
    cmd = [
        "gh",
        "pr",
        "list",
        "--base",
        target_branch,
        "--state",
        "merged",
        "--limit",
        "1000",
        "--json",
        "number,createdAt,mergedAt",
    ]
    output = subprocess.check_output(cmd, text=True)
    return json.loads(output)


def is_model_changed_in_pr(pr_number):
    # Get the list of files changed in the PR.
    # We use the GitHub API (paginated) rather than `gh pr diff --name-only`
    # because `gh pr diff` fetches the full patch and fails on PRs with very
    # large diffs (returning a non-zero exit / "stream CANCEL"). The API's
    # pulls/{number}/files endpoint just lists filenames and paginates, so it
    # handles large PRs gracefully.
    cmd = [
        "gh",
        "api",
        f"repos/:owner/:repo/pulls/{pr_number}/files",
        "--paginate",
        "--jq",
        ".[].filename",
    ]
    try:
        output = subprocess.check_output(cmd, text=True)
    except subprocess.CalledProcessError as e:
        # Don't let one problematic PR crash the whole run; warn and treat it
        # as "not changed" so it is skipped.
        warnings.warn(
            f"Could not fetch changed files for PR #{pr_number} "
            f"(exit {e.returncode}); skipping it."
        )
        return False
    changed_files = output.splitlines()
    # Check if 'model.xml' is in the list of changed files
    return "model.xml" in changed_files


if __name__ == "__main__":
    run_tests_on_prs(
        pr_start=PR_START,
        pr_end=PR_END,
        flux_limit=DEFAULT_FLUX_LIMIT,
        biomass_rxn_id=BIOMASS_RXN_ID,
        force_rerun=FORCE_RERUN,
    )
    print("You can plot the results using plot_match_over_time.py")

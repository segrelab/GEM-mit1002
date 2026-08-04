import json
import os
import subprocess
import warnings

import cobra
import pandas as pd
from gem_utilities import media

# ---- CONFIG: edit these, then run `python run_tests_on_prs.py` ----
# Range of merged PRs to evaluate. PR_END = None goes up to the latest PR.
# For the initial model construction to growing on everything use 89 - 212
# To highlight the acetate/leucine/isoleucine fixes use 273 - 293
# For initial to after fixing the TCA cycle fluxes use 89 - 344
PR_START = 89
PR_END = None
# Flux magnitude above which a reaction is considered problematically unbounded
UNBOUNDED_FLUX_LIMIT = 100
# Reaction ID for the biomass reaction
BIOMASS_RXN_ID = "bio1_biomass"
# False = incremental (only new PRs, plus any that previously errored).
# True = re-run every PR in the range, overwriting existing results
# (e.g. when the test itself has changed).
FORCE_RERUN = False

# FILE PATHS
FILE_PATH = os.path.dirname(os.path.abspath(__file__))
REPO_PATH = os.path.dirname(FILE_PATH)

import sys

sys.path.insert(0, str(REPO_PATH))

from tools.media import MEDIA  # noqa: E402
DATA_DIR = os.path.join(REPO_PATH, "data")

# Load the media definitions
media_definitions = MEDIA

# Load the TSV of the growth phenotypes
growth_phenotypes = pd.read_csv(
    os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"),
    sep="\t",
    converters={"met_id": lambda x: x.split(",")},
)


def run_tests_on_prs(
    pr_start=89,
    pr_end=None,
    unbounded_flux_limit: int = 999,
    biomass_rxn_id="bio1_biomass",
    force_rerun: bool = False,
    prs_to_skip=[],
):
    """Run the growth tests on merged PRs and save a summary CSV.

    By default this is incremental: PRs that already have a (non-ERROR) row in
    growth_match_summary.csv are skipped, so only new PRs (and any that
    previously errored) are evaluated. Existing results are preserved and the
    new ones are merged in.

    Set ``force_rerun=True`` to re-evaluate every PR in the range, overwriting
    the existing rows for those PRs (useful when the test itself has changed).
    Results for PRs outside the requested range are always preserved.
    """
    # Load any existing results so we can run incrementally. Keyed by PR number.
    summary_file = os.path.join(FILE_PATH, "growth_match_summary.csv")
    existing_results = {}
    if os.path.exists(summary_file):
        existing_df = pd.read_csv(summary_file)
        existing_results = {
            int(row["PR Number"]): row.to_dict() for _, row in existing_df.iterrows()
        }
    # PRs we already have a successful (non-ERROR) result for. Unless
    # force_rerun is set, these are skipped. PRs whose stored result is ERROR
    # are always re-attempted.
    already_done = {
        pr_number
        for pr_number, row in existing_results.items()
        if str(row.get("Matches")) != "ERROR"
    }

    # Get a list of PRs merged into the dev branch (could change to main, by
    # providing target_branch="main")
    # By only looking at dev we miss a few PRs that were merged into main,
    # (89-117) but this is easier to automate and still captures the majority
    # of changes to the model
    all_pr_entries = get_prs_by_target()
    # Filter PRs for a specific range of interest
    # For the initial model construction to growing on everything use PRs 89 - 212
    # To highlight the acetate/leucine/isolecuine fixes use PRs 273-293
    # For initial to after fixing the TCA cycle fluxes use PRs 89-344
    # If no upper limit to the range is given, go to the latest version
    if pr_end is None:
        pr_end = max(pr_entry["number"] for pr_entry in all_pr_entries)
    pr_entries_to_check = [
        pr_entry
        for pr_entry in all_pr_entries
        if pr_start <= pr_entry["number"] <= pr_end
        and pr_entry["number"] not in prs_to_skip
        # Incremental: skip PRs we already have results for (unless forcing)
        and (force_rerun or pr_entry["number"] not in already_done)
    ]
    # Filter PRs to those that changed the model.xml file
    pull_requests = [
        pr_entry
        for pr_entry in pr_entries_to_check
        if is_model_changed_in_pr(pr_entry["number"])
    ]

    if not pull_requests:
        print("No new PRs to evaluate; existing results are up to date.")
    else:
        print(f"Evaluating {len(pull_requests)} PR(s)...")

    # Start from the existing results and update them with anything we (re)run
    results_by_pr = dict(existing_results)

    for pr in pull_requests:
        print(f"\n--- Evaluating PR #{pr['number']} ---")

        try:
            # Use gh to get the content of model.xml for a specific PR
            subprocess.run(
                [
                    "gh",
                    "api",
                    f"repos/:owner/:repo/contents/model.xml?ref=pull/{pr['number']}/head",
                    "-H",
                    "Accept: application/vnd.github.v3.raw",
                ],
                stdout=open("temp_model.xml", "w"),
            )

            # Run test_growth
            results = run_test_growth(
                unbounded_flux_limit=unbounded_flux_limit, biomass_rxn_id=biomass_rxn_id
            )

            # Store results (overwrites any existing row for this PR)
            results_by_pr[pr["number"]] = {
                "PR Number": pr["number"],
                "Date Opened": pr["createdAt"],
                "Date Merged": pr["mergedAt"],
                "Reactions": results.get("num_reactions"),
                "Metabolites": results.get("num_metabolites"),
                "Genes": results.get("num_genes"),
                "Matches": results.get("matches"),
                "Total": results.get("total"),
                "% Match": round(
                    100 * results.get("matches", 0) / max(results.get("total", 1), 1),
                    2,
                ),
                "Unbounded Flux Reactions": results.get("num_unbounded_rxns"),
            }

        except Exception as e:
            print(f"Error with PR #{pr['number']}: {e}")
            results_by_pr[pr["number"]] = {
                "PR Number": pr["number"],
                "Date Opened": pr["createdAt"],
                "Date Merged": pr["mergedAt"],
                "Reactions": "ERROR",
                "Metabolites": "ERROR",
                "Genes": "ERROR",
                "Matches": "ERROR",
                "Total": "ERROR",
                "% Match": "ERROR",
                "Unbounded Flux Reactions": "ERROR",
            }

    # Delete the temporary model file
    if os.path.exists("temp_model.xml"):
        os.remove("temp_model.xml")

    # Write results to CSV ONCE, after all PRs are processed. Sort by PR number
    # so newly added PRs land in order alongside the existing ones.
    results_list = [results_by_pr[num] for num in sorted(results_by_pr)]
    df = pd.DataFrame(results_list)
    df.to_csv(summary_file, index=False)


def run_test_growth(unbounded_flux_limit: int = 999, biomass_rxn_id="bio1_biomass"):
    """
    On the current branch, this function re-runs the growth with pFBA, while
    holding the amount carbon constant across sources, and saves how many
    media conditions support growth, and how many reactions across all
    simulations have a flux above an arbitrary threshold.

    Parameters
    ----------
    unbounded_flux_limit : int, optional
        The value at which a flux is considered problematically large, by
        default 999
    biomass_rxn_id : str, optional
        Reaction ID for the biomass reaction, by default "bio1_biomass"

    Returns
    -------
    dict
       Dictionary of the results for the model, including the number of
       reactions, metabolites, and genes, the number of matches between
       predicted and expected growth phenotypes, the total number of phenotypes
       tested, and the number of unique reactions with a flux above the
       unbounded_flux_limit across all simulations.
    """
    model = cobra.io.read_sbml_model(os.path.join(REPO_PATH, "temp_model.xml"))

    # Count model components
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)

    # Start counters
    matches = 0
    total = 0
    unique_rxns_with_unbounded_flux = set()

    for _, row in growth_phenotypes.iterrows():
        # Get the expected growth phenotype
        expected_growth = row["growth"]

        # If the expected growth is not Yes or No (e.g. "Unsure"), skip the row
        if expected_growth not in ["Yes", "No"]:
            warnings.warn(
                f"Expected growth phenotype '{expected_growth}' is not valid. Skipping row."
            )
            continue

        # Set the minimal media and exchange reactions
        minimal_media = media_definitions[row["minimal_media"]].copy()

        # Add the metabolite(s) specified in the row to the media
        for met_id in row["met_id"]:
            # Check that there is an exchange reaction for the metabolite in the model
            if "EX_" + met_id + "_e0" not in [r.id for r in model.reactions]:
                warnings.warn(f"Model does not have an exchange reaction for {met_id}.")
                continue
            # Get the metabolite object
            met = model.metabolites.get_by_id(met_id + "_e0")
            # Get the number of carbon atoms in the metabolite
            n_carbons = met.elements.get("C", 0)
            # If the metabolite has carbons, set the lower bound to be 60/n_carbons (equivalent to 10 for glucose)
            # If the metabolite does not have carbon, set an unlimited amount
            if n_carbons == 0:
                minimal_media["EX_" + met_id + "_e0"] = 1000.0
            else:
                minimal_media["EX_" + met_id + "_e0"] = 60 / n_carbons
        # Set the media
        model.medium = media.clean_media(model, minimal_media)
        # Run pFBA on the model.
        # Since the model has a non-zero maintenance requirement (a forced
        # lower bound on ATP hydrolysis), a medium that cannot supply that
        # energy makes the LP infeasible and pfba raises. Biologically that
        # just means the organism can't sustain itself, i.e. no growth, so we
        # catch it and treat it as such rather than letting it error out the
        # whole PR.
        try:
            sol = cobra.flux_analysis.pfba(model)
        except (cobra.exceptions.Infeasible, cobra.exceptions.OptimizationError):
            sol = None
        # Check if the model grows
        if sol is not None and sol.fluxes[biomass_rxn_id] > 1e-3:
            # If it does, set to Yes
            pred_growth = "Yes"
            # Get the number of reactions in the solution with a flux above the unbounded_flux_limit
            rxns_with_unbounded_flux = [
                r for r, flux in sol.fluxes.items() if abs(flux) > unbounded_flux_limit
            ]
            # Add the unique reactions with unbounded fluxes to the set
            unique_rxns_with_unbounded_flux.update(rxns_with_unbounded_flux)
        else:
            # If it doesn't (or the problem was infeasible), set to No
            pred_growth = "No"

        if pred_growth == expected_growth:
            matches += 1
        total += 1

    return {
        "num_reactions": num_reactions,
        "num_metabolites": num_metabolites,
        "num_genes": num_genes,
        "matches": matches,
        "total": total,
        "num_unbounded_rxns": len(unique_rxns_with_unbounded_flux),
    }


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
        unbounded_flux_limit=UNBOUNDED_FLUX_LIMIT,
        biomass_rxn_id=BIOMASS_RXN_ID,
        force_rerun=FORCE_RERUN,
    )
    print("Test results saved to growth_match_summary.csv")
    print("You can plot the results using the provided plotting script.")

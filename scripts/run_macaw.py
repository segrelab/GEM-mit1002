"""Run the MACAW test suite against the model.

Writes ``macaw_results.csv`` and ``macaw_edge_list.csv`` to ``scripts/results/``
alongside the other generated reports.

This is *not* run in CI, because the dilution test takes too long. Run it by hand
from anywhere:

    python scripts/run_macaw.py
"""

import os

import cobra
from macaw.main import run_all_tests

# import py4cytoscape as p4c

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(REPO_ROOT, "model.xml")
RESULTS_DIR = os.path.join(REPO_ROOT, "scripts", "results")


# Run the MACAW pipeline
def run_macaw(model):
    # TODO: Figure out how to run the function without having to wrap it
    # Run all tests
    (test_results, edge_list) = run_all_tests(model)

    return test_results, edge_list


if __name__ == "__main__":
    # Load the model
    model = cobra.io.read_sbml_model(MODEL_PATH)

    # Run the MACAW pipeline
    (test_results, edge_list) = run_all_tests(model)

    # Save the results
    os.makedirs(RESULTS_DIR, exist_ok=True)
    test_results.to_csv(os.path.join(RESULTS_DIR, "macaw_results.csv"))
    edge_list.to_csv(os.path.join(RESULTS_DIR, "macaw_edge_list.csv"))

    # TODO: Visualize the network results

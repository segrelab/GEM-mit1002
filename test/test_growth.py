import os
import unittest
import warnings

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from gem_utilities import biomass, media
from tools.media import MEDIA

# Set path to the `test_files` directory
TESTFILE_DIR = os.path.join(os.path.dirname(__file__), "test_files")
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

# Set path to the results directory
RESULTS_DIR = os.path.join(TESTFILE_DIR, "test_results")

# Load the media definitions
media_definitions = MEDIA
minimal_media = media_definitions["minimal"]


class TestGrowthPhenotypes(unittest.TestCase):
    # Check that there is no growth on a media with no carbon sources
    def test_growth_w_0_C(self):
        # Load the model with cobrapy
        model = cobra.io.read_sbml_model("model.xml")

        # Set the media so that there are no carbon sources
        model.medium = media.clean_media(model, minimal_media)

        # Run the model
        sol = model.optimize()

        # Check that no biomass is produced
        self.assertEqual(sol.objective_value, 0)

    # Test that there is growth or no growth as expected on different media
    def test_expected_growth_phenotypes(self):
        # Load the TSV of the growth phenotypes
        growth_phenotypes = pd.read_csv(
            os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"),
            sep="\t",
            converters={"met_id": lambda x: x.split(",")},
        )

        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Loop through the growth phenotpes, and add the carbon source to the
        # minimal media, run FBA and check if the model grows
        ex_rxn_present = []
        pred_growth = []
        for index, row in growth_phenotypes.iterrows():
            minimal_media = media_definitions[row["minimal_media"]].copy()
            # Check if the model has an exchange reaction for the metabolite(s)
            if all(
                "EX_" + met_id + "_e0" in [r.id for r in model.reactions]
                for met_id in row["met_id"]
            ):
                # If it does, add the exchange reaction to the minimal media used
                for met_id in row["met_id"]:
                    minimal_media["EX_" + met_id + "_e0"] = 1000.0
                # Mark the exchange reaction as present
                ex_rxn_present.append("Yes")
            else:
                # Mark the exchange reaction as not present
                ex_rxn_present.append("No")
            # Set the media
            model.medium = media.clean_media(model, minimal_media)
            # Run the model
            sol = model.optimize()
            # Check if the model grows
            if sol.objective_value > 1e-3:
                # If it does, add 'Y' to the list
                pred_growth.append("Yes")
            else:
                # If it doesn't, add 'N' to the list
                pred_growth.append("No")

        # Add the lists as new columns in the dataframe
        growth_phenotypes["all_ex_rxn_present"] = ex_rxn_present
        growth_phenotypes["pred_growth"] = pred_growth

        # Filter for only the phenotypes that are explicitly "Yes" or "No"
        # Since these are the only ones we can test
        testable_phenotypes = growth_phenotypes[
            growth_phenotypes["growth"].isin(["Yes", "No"])
        ].copy()

        # Find all rows where the experimental and predicted growth do not match
        mismatches = testable_phenotypes[
            testable_phenotypes["growth"] != testable_phenotypes["pred_growth"]
        ]

        # If the mismatches DataFrame is not empty, the test fails.
        if not mismatches.empty:
            # Format a detailed error message
            error_message = [
                "\nPredicted growth phenotypes do not match experimental data for the following conditions:"
            ]
            # Add a header for the output table
            header = f"{'Media':<35} | {'Carbon Source':<25} | {'Experimental':<12} | {'Predicted':<12} | {'Ex Rxn Present':<10}"
            error_message.append(header)
            error_message.append("-" * len(header))

            # Add a row for each mismatch
            for index, row in mismatches.iterrows():
                # Note: Assuming you have a 'c_source' column as in your original file
                line = f"{row['minimal_media']:<35} | {row['c_source']:<25} | {row['growth']:<12} | {row['pred_growth']:<12} | {row['all_ex_rxn_present']:<10}"
                error_message.append(line)

            # Combine the list into a single string and fail the test
            self.fail("\n".join(error_message))


if __name__ == "__main__":
    unittest.main()

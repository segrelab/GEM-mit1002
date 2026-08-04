import os
import warnings

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(FILE_DIR)

import sys

sys.path.insert(0, str(REPO_DIR))

from tools.media import MEDIA  # noqa: E402
DATA_DIR = os.path.join(REPO_DIR, "data")
# Define the directory containing the files
gapfill_dir = os.path.join(REPO_DIR, "modelseedpy_gapfill_per_biomass_cmpt")

OUT_DIR = os.path.join(FILE_DIR, "component_producibility_results")
# If the output directory doesn't exist, create it
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)


def main():
    # Make a dctionary of the model IDs and the file paths
    model_files = {
        # "Base Model (ModelSEEDpy)": os.path.join(REPO_DIR, "modelseedpy_model_01.xml"),
        # "Glucose Gapfilled (ModelSEEDpy)": os.path.join(
        #     REPO_DIR, "modelseedpy_model_04.xml"
        # ),
        "KBase Model, ModelSEEDpy Gapfilled": os.path.join(
            REPO_DIR, "2025-01-08_Scott_draft-model-from-KBase-MSP-gapfilled.xml"
        ),
    }

    # # Load the base model into COBRA
    # # Needed to get the metabolite names
    # base_model = cobra.io.read_sbml_model(model_files["Base Model (ModelSEEDpy)"])

    # # Add all of the models that were gpafilled (on mbm + glucose) for a single biomass component to the dictionart
    # # Iterate through files and add to dictionary
    # for filename in os.listdir(gapfill_dir):
    #     if filename.startswith("gapfilled_cpd") and filename.endswith(".xml"):
    #         # Extract cpdXXXXX_c0 from the filename (remove "gapfilled_" and ".xml")
    #         compound_id = filename[10:-4]
    #         # Get the human-friendly name of the metabolite
    #         met_name = base_model.metabolites.get_by_id(compound_id).name
    #         key = f"Gapfilled for {met_name}"
    #         model_files[key] = os.path.join(gapfill_dir, filename)

    # Load the models and update their IDs (used in the result filenames) and store in a list for easy access
    models = []
    for model_id, model_file in model_files.items():
        # Load the model
        model = cobra.io.read_sbml_model(model_file)
        model.id = model_id
        models.append(model)

    # Load the growth pheonotype results
    growth_phenotypes = pd.read_csv(
        os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"),
        sep="\t",
        header=0,
    )

    # Filter the growth phenotypes to only include the carbon sources that it can grow on
    growth_phenotypes = growth_phenotypes[growth_phenotypes["growth"] == "Yes"]

    # Load the media definitions
    media_definitions = MEDIA

    # Run my function on each of the models
    for model in models:
        test_model(
            model, growth_phenotypes, media_definitions, biomass_rxn="bio1_biomass"
        )


# Helper function for setting the media regardless if the exchange reaction is
# present in the model
# TODO: Move this to a helper file, and remove from the test growth file too
def clean_media(model, media):
    """clean_media
    Removes exchange reactions from the media that are not present in the model

    Parameters
    ----------
    model : cobra.Model
        The model to set the media for.
    media : dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction.

    Returns
    -------
    dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction
    """
    # Make an empty dictionary for the media
    clean_medium = {}
    # Loop through the media and set the exchange reactions that are present
    for ex_rxn, lb in media.items():
        if ex_rxn in [r.id for r in model.reactions]:
            clean_medium[ex_rxn] = lb
        else:
            warnings.warn(
                "Model does not have the exchange reaction "
                + ex_rxn
                + ", so it was not set in the media."
            )

    # Return the clean medium
    return clean_medium


# Broke out the function for actually testing the component producibility on a given medium
def test_medium(medium_dict, biomass_compounds, model):
    """Runs FBA for each biomass_compound as objective under `medium_dict`.
    Returns a dict {compound_id: True/False}."""
    results = {}
    # Set the medium
    model.medium = clean_media(model, medium_dict)

    # Loop over each biomass compound
    for cpd_id in biomass_compounds:
        # Get the human-friendly name of the metabolite
        cpd_name = model.metabolites.get_by_id(cpd_id).name

        # Set the objective to the sink/demand for cpd_id
        # e.g., "SK_cpd_id"
        model.objective = {model.reactions.get_by_id("SK_" + cpd_id): 1}

        # Optimize
        sol = model.optimize()

        # Check feasibility
        if (sol.status == "optimal") and (sol.objective_value > 1e-6):
            results[cpd_name] = True
        else:
            results[cpd_name] = False
    return results


def test_model(model, growth_phenotypes, media_definitions, biomass_rxn="bio1_biomass"):
    """
    Add exchange reactions for all metabolites in the model and test the producibility of the biomass components on the given media

    Args:
    model (cobra.Model): The model to test
    biomass_rxn (str): The ID of the biomass reaction in the model. Default is "bio1_biomass" (which is used in KBase models).
    """
    # Check that the growth phenotypes dataframe has the expected columns
    expected_columns = ["minimal_media", "c_source", "met_id", "growth"]
    if not all(col in growth_phenotypes.columns for col in expected_columns):
        raise ValueError(
            "growth_phenotypes dataframe must have columns "
            + ", ".join(expected_columns)
        )

    # Add sink reactions for all metabolites, but set the lower bound to 0
    # because by default the sink reactions are reversible, and so can be
    # used to import metabolites that are not in the media
    for metabolite in model.metabolites:
        # Check if there is already a sink reaction for this metabolite
        if "SK_" + metabolite.id not in [r.id for r in model.reactions]:
            model.add_boundary(metabolite, type="sink", lb=0)

    # Get the biomass composition from the model
    biomass_rxn = model.reactions.get_by_id(biomass_rxn)
    biomass_compounds = [
        met.id for met in biomass_rxn.metabolites if biomass_rxn.metabolites[met] < 0
    ]

    # Make a dictionary to store the producibility results
    biomass_producibility = {}

    # Add negative controls
    # TODO: Create the list of controls programmatically from the growth phenotypes data
    # Create some custom dicts for the "control" conditions:
    # 1) An empty medium
    empty_media = {}

    # 2) mbm minus carbon sources (assuming your minimal media has some "EX_CARBON" keys)
    #    This just filters out any exchange that might be for the primary carbon.
    mbm_no_carbon = {
        ex: lb
        for ex, lb in media_definitions["mbm_media"].items()
        if not ex.startswith("EX_") or "glc" not in ex.lower()
    }

    # 3) l1 minus carbon sources
    l1_no_carbon = {
        ex: lb
        for ex, lb in media_definitions["l1_media"].items()
        if not ex.startswith("EX_") or "glc" not in ex.lower()
    }

    # Combine the negative controls into a dictionary
    controls = {"empty": empty_media, "mbm_noC": mbm_no_carbon, "l1_noC": l1_no_carbon}

    # Test the negative controls
    for ctrl_name, ctrl_medium in controls.items():
        biomass_producibility[ctrl_name] = test_medium(
            ctrl_medium, biomass_compounds, model
        )

    # Check the producibility of the biomass componets on the different carbon sources
    for index, row in growth_phenotypes.iterrows():
        # Make an ID for the results that is combination of the minimal media name and the carbon source
        c_source = row["minimal_media"] + "_" + row["c_source"]
        # Make a dictionary to store the results for just this carbon source
        biomass_producibility[c_source] = {}
        # Set the model media to match the experimental media
        medium = media_definitions[row["minimal_media"]].copy()
        medium["EX_" + row["met_id"] + "_e0"] = (
            1000.0  # FIXME: I should set this to a consistent, lower value
        )
        # Test it
        biomass_producibility[c_source] = test_medium(medium, biomass_compounds, model)

    # Make a dataframe of the producibility results and save it to a CSV file
    df = pd.DataFrame.from_dict(biomass_producibility)
    # Save the dataframe to a CSV file and make the file name specific the the model.id
    df.to_csv(os.path.join(OUT_DIR, model.id + "_biomass_producibility.csv"))

    # Plot the producibility results
    plot_prodcubility(model, df)


# Function for plotting the producibility results
def plot_prodcubility(model, df):
    # Convert the producibility results to a boolean
    df_binary = df.replace({True: 1, False: 0})

    # Define my own color palette
    my_palette = sns.color_palette(["red", "green"])

    # Plot the producibility results as a heatmap
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(
        df_binary,
        annot=False,
        cmap=my_palette,
        cbar_kws={"ticks": [0.25, 0.75], "label": "Producibile"},
    )

    # Now fix the colorbar tick labels
    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(["No", "Yes"])

    # Add white lines to separate the different cells
    for i in range(df_binary.shape[0]):
        plt.axhline(i, color="white", linewidth=0.5)
    for i in range(df_binary.shape[1]):
        plt.axvline(i, color="white", linewidth=1)

    # Make sure all y ticks/component names are shown
    plt.yticks(
        [x + 0.5 for x in range(df_binary.shape[0])], df_binary.index, fontsize=5
    )

    plt.xlabel("Biomass composition")
    plt.ylabel("Biomass compound")

    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(OUT_DIR, model.id + "_biomass_producibility_heatmap.png"))


if __name__ == "__main__":
    main()

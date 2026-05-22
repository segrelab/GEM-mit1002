import os
import pickle

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from gem_utilities import biomass, media

# Define paths relative to the script or project root
# It's better practice to define a project root
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
TESTFILE_DIR = os.path.join(PROJECT_ROOT, "test", "test_files")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "scripts", "results")

# Load the media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)

# Define a dictionary of human-friendly versions of the media names
# The key is the name in the media_definitions dictionary
# The value is the human-friendly name to use the table
media_names = {
    "l1": "L1",
    "mbm": "Minimal Basal Medium (Moran Lab)",
    "promm_no_c": "ProMM",
    "marine_broth_wo_yeast_and_peptone": "Marine Broth",
    "marine_broth_wo_yeast_and_peptone_no_n": "Marine Broth (No Nitrogen)",
    "swm": "Seawater Medium",
}


def generate_growth_phenotype_report(model: cobra.Model):
    # Load the TSV of the growth phenotypes
    growth_phenotypes = pd.read_csv(
        os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"),
        sep="\t",
        converters={"met_id": lambda x: x.split(",")},
    )

    # Loop through the growth phenotpes, and add the carbon source to the
    # minimal media, run FBA and check if the model grows
    ex_rxn_present = []
    fba_growth_rate = []
    pred_growth = []
    for index, row in growth_phenotypes.iterrows():
        minimal_media = media_definitions[row["minimal_media"]].copy()
        # Check if the model has an exchange reaction for the metabolite
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
        # Save the growth rate
        fba_growth_rate.append(sol.objective_value)
        # Save the
        if sol.objective_value > 1e-3:
            # If it does, add 'Y' to the list
            pred_growth.append("Yes")
        else:
            # If it doesn't, add 'N' to the list
            pred_growth.append("No")

    # Add the lists as new columns in the dataframe
    growth_phenotypes["all_ex_rxn_present"] = ex_rxn_present
    growth_phenotypes["pred_growth"] = pred_growth
    growth_phenotypes["fba_growth_rate"] = fba_growth_rate

    # Save the dataframe as a TSV
    growth_phenotypes.to_csv(
        os.path.join(RESULTS_DIR, "known_growth_phenotypes_w_pred.tsv"),
        sep="\t",
        index=False,
    )

    # Plot a categorical heatmap of the growth phenotypes, where the rows
    # are the metabolites and the columns are the experimental and predicted
    # growth phenotypes. Show growth as blue and no growth as orange, and
    # unsure as gray
    # First, make a new dataframe with the metabolites as the rows and the
    # experimental and predicted growth phenotypes as the columns
    # Combine the values of "minimal_media" and "c_source" into one column
    growth_phenotypes["c_source"] = (
        growth_phenotypes["minimal_media"] + " " + growth_phenotypes["c_source"]
    )
    # And set it as the index
    growth_phenotypes = growth_phenotypes.set_index("c_source")

    # Make a dictionary for the phenotypes to numbers
    value_to_int = {"Unsure": 0, "No": 1, "Yes": 2}
    n = len(value_to_int)

    # Create an annotation data frame for the text labels on the heatmap
    annotation_key = {"No": "No Exchange", "Yes": ""}
    annot_df = (
        growth_phenotypes["all_ex_rxn_present"].replace(annotation_key).to_frame()
    )
    annot_df.rename(columns={"all_ex_rxn_present": "FBA"}, inplace=True)
    annot_df["Experimental"] = ""
    # Sort the columns to match the order of the heatmap
    annot_df = annot_df[["Experimental", "FBA"]]

    # Subset the other columns, to have just the growth and predicted growth
    growth_phenotypes = growth_phenotypes[["growth", "pred_growth"]]

    # Rename the columns and the index to be longer/more descriptive
    growth_phenotypes.index.name = "Media/Carbon Source"
    growth_phenotypes = growth_phenotypes.rename(
        columns={
            "growth": "Experimental",
            "pred_growth": "FBA",
        }
    )

    # Make a colormap of specified colors (in numerical order for the phenotypes)
    # cmap = ['gray', '#F18F01', '#399E5A'] # Gray, orange, green
    cmap = ["#5E5E5E", "#FF7D0A", "#024064"]  # C-CoMP gray, orange, and dark blue

    # Dynamically set the figure height based on the number of rows
    fig_height = max(10, len(growth_phenotypes) * 0.4)  # 0.4 inches per row
    # Plot the heatmap
    # Use constrained_layout to prevent cutting off y-axis/colorbar labels
    fig, ax = plt.subplots(
        figsize=(8, fig_height),
        constrained_layout=True,
    )
    sns.heatmap(
        growth_phenotypes.replace(value_to_int),
        cmap=cmap,
        linewidths=4,
        linecolor="white",
        annot=annot_df,
        fmt="",
        annot_kws={"fontsize": 8},  # Smaller font size for annotation
        ax=ax,
    )

    # Modify colorbar:
    colorbar = ax.collections[0].colorbar
    r = colorbar.vmax - colorbar.vmin
    colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
    colorbar.set_ticklabels(list(value_to_int.keys()))

    # Move the x-axis labels to the top
    plt.tick_params(
        axis="both",
        which="major",
        labelsize=10,
        labelbottom=False,
        bottom=False,
        top=True,
        labeltop=True,
    )

    # Make sure that every y-tick is shown
    ax.set_yticks([i + 0.5 for i in range(len(growth_phenotypes))])
    ax.set_yticklabels(growth_phenotypes.index, rotation=0)

    # Save the figure
    plt.savefig(os.path.join(RESULTS_DIR, "exp_vs_pred_growth_phenotypes.png"))


def generate_biomass_producibility_report(model: cobra.Model):
    # Load the TSV of the growth phenotypes
    growth_phenotypes = pd.read_csv(
        os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"),
        sep="\t",
        converters={"met_id": lambda x: x.split(",")},
    )

    # Filter the growth phenotypes to only include the carbon sources that it can grow on
    growth_phenotypes = growth_phenotypes[growth_phenotypes["growth"] == "Yes"]

    # Run the biomass producibility function on each of the models
    sink_options = [False, True]
    for option in sink_options:
        biomass.check_biomass_producibility(
            model,
            growth_phenotypes,
            media_definitions,
            sinks_for_all=option,
            out_dir=RESULTS_DIR,
        )


if __name__ == "__main__":
    # Ensure the results directory exists
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(PROJECT_ROOT, "model.xml"))

    # Generate the reports
    generate_growth_phenotype_report(model)
    generate_biomass_producibility_report(model)

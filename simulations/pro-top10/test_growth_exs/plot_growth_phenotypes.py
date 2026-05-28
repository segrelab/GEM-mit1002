import pickle
from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from gem2cue import (
    utils,  # Import the working version (works with the med4-hot1a3 conda env)
)

FILE_DIR = Path(__file__).resolve().parent
TOP_10_DIR = FILE_DIR.parent
REPO_DIR = FILE_DIR.parents[2]

# Set path to the `test_files` directory
TESTFILE_DIR = REPO_DIR / "test" / "test_files"

# Load the media definitions
with open(TESTFILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
    media_definitions = pickle.load(f)
minimal_media = media_definitions["minimal"]

# Load the model with cobrapy
model_orig = cobra.io.read_sbml_model(REPO_DIR / "model.xml")
c_ex_rxns = utils.get_c_ex_rxns(model_orig)

# Load the top 10 metabolite file
top_10_metabolites = pd.read_csv(
    TOP_10_DIR / "data" / "top10_exometabolites.csv", sep=","
)

# Load the known growth phenotype file
known_growth_phenotypes = pd.read_csv(
    TESTFILE_DIR / "known_growth_phenotypes.tsv", sep="\t"
)

# Subset the top 10 metabolites to only include the ones that are from PRO
top_10_metabolites = top_10_metabolites[
    top_10_metabolites["organism"] == "Prochlorococcus marinus"
].copy()

# For every metabolite in the top 10 metabolites, check if it is in the
# known growth phenotypes, if it is add add the value from that file to the
# new column "exp_growth" in the top 10 metabolites dataframe, then run FBA
# and check if the model grows on that metabolite, and add the value to the
# new column "pred_growth" in the top 10 metabolites dataframe
for index, row in top_10_metabolites.iterrows():
    # Check if the metabolite is in the known growth phenotypes
    if row["met_id"] in known_growth_phenotypes["met_id"].values:
        # Check that all the values in the known growth file are the same
        exp_growth = known_growth_phenotypes[
            known_growth_phenotypes["met_id"] == row["met_id"]
        ]["growth"].values
        if len(set(exp_growth)) > 1:
            # As long as there is one "Yes" in the list, set the value to "Yes"
            if "Yes" in exp_growth:
                top_10_metabolites.at[index, "exp_growth"] = "Yes"
            # If there is no "Yes" in the list, throw an error
            else:
                raise ValueError(
                    f"Multiple growth phenotypes found for metabolite {row['met_id']}: {exp_growth}"
                )
        else:
            # If there is only one value, set the value to that
            top_10_metabolites.at[index, "exp_growth"] = exp_growth[0]
    else:
        # If it isn't, set the value to "Unknown"
        top_10_metabolites.at[index, "exp_growth"] = "Unknown"

    # Make a copy of the model to work with
    model = model_orig.copy()

    # Set the media so that there are no carbon sources
    model.medium = minimal_media

    # Check if the model has an exchange reaction for the metabolite
    if "EX_" + row["met_id"] + "_e0" in [r.id for r in model.reactions]:
        # If it does, add the exchange reaction to the minimal media used
        medium = model.medium.copy()
        medium["EX_" + row["met_id"] + "_e0"] = 1000.0
        model.medium = medium
        # Say "Yes" to the exchange reaction present
        top_10_metabolites.at[index, "ex_rxn_present"] = "Yes"
    # If not, save that information in a new column
    else:
        top_10_metabolites.at[index, "ex_rxn_present"] = "No"

    # Run the model
    sol = model.optimize()

    # Check if the model grows
    if sol.objective_value > 1e-3:
        # If it does, add 'Y' to the list
        top_10_metabolites.at[index, "pred_growth"] = "Yes"

        # Extract the carbon fates for the solution (not normalized)
        c_fates = utils.extract_c_fates_from_solution(
            sol, c_ex_rxns, co2_ex_rxn="EX_cpd00011_e0", norm=False
        )
        uptake = c_fates[0]
        co2 = c_fates[1]
        organic_c = c_fates[2]
        biomass = c_fates[3]

        # Calculate CUE/GGE/BGE
        cue = 1 - (co2 / uptake)
        gge = biomass / uptake
        bge = biomass / (biomass + co2)

        # Add the results to the dataframe
        top_10_metabolites.at[index, "cue"] = cue
        top_10_metabolites.at[index, "gge"] = gge
        top_10_metabolites.at[index, "bge"] = bge
    else:
        # If it doesn't, add 'N' to the list
        top_10_metabolites.at[index, "pred_growth"] = "No"
        # Set the results to NaN
        top_10_metabolites.at[index, "cue"] = None
        top_10_metabolites.at[index, "gge"] = None
        top_10_metabolites.at[index, "bge"] = None

    # Check for growth with free exchnage and tranport (i.e. add a sink reaction,
    # and check if the model grows)
    # Make a new copy of the model to work with
    model = model_orig.copy()
    # Set the media so that there are no carbon sources
    model.medium = minimal_media
    # Add a sink reaction for the metabolite
    # TODO: Change the lower bound to restrict the flux
    model.add_boundary(
        metabolite=model.metabolites.get_by_id(row["met_id"] + "_c0"),
        type="sink",
    )
    # Run the model
    sol = model.optimize()
    # Check if the model grows
    if sol.objective_value > 1e-3:
        # If it does, add 'Y' to the list
        top_10_metabolites.at[index, "pred_growth_free_transport"] = "Yes"
    else:
        # If it doesn't, add 'N' to the list
        top_10_metabolites.at[index, "pred_growth_free_transport"] = "No"

# Plot a categorical heatmap of the growth phenotypes, where the rows
# are the metabolites and the columns are the experimental and predicted
# growth phenotypes. Show growth as blue and no growth as red, and
# unsure as gray
# And set it as the index
growth_phenotypes = top_10_metabolites.set_index("metabolite")

# Make a dictionary for the phenotypes to numbers
value_to_int = {"Unknown": 0, "No": 1, "Yes": 2}
n = len(value_to_int)

# Create an annotation data frame for the text labels on the heatmap
annotation_key = {"No": "No Exchange", "Yes": ""}
annot_df = growth_phenotypes["ex_rxn_present"].replace(annotation_key).to_frame()
annot_df.rename(columns={"ex_rxn_present": "FBA (Defined EXs Only)"}, inplace=True)
annot_df["Experimental"] = ""
annot_df["FBA (Free Transport)"] = ""
# Sort the columns to match the order of the heatmap
annot_df = annot_df[["Experimental", "FBA (Defined EXs Only)", "FBA (Free Transport)"]]

# Subset the other columns, to have just the growth and the two types of
# predicted growth (with the model as is, and with free transport)
growth_phenotypes_to_plot = growth_phenotypes[
    ["exp_growth", "pred_growth", "pred_growth_free_transport"]
]

# Rename the columns and the index to be longer/more descriptive
growth_phenotypes_to_plot.index.name = "Media/Carbon Source"
growth_phenotypes_to_plot = growth_phenotypes_to_plot.rename(
    columns={
        "exp_growth": "Experimental",
        "pred_growth": "FBA (Defined EXs Only)",
        "pred_growth_free_transport": "FBA (Free Transport)",
    }
)

# Make a colormap of specified colors (in numerical order for the phenotypes)
cmap = ["#5E5E5E", "#AC333C", "#193F61"]  # C-CoMP gray, red, and dark blue

# Plot the heatmap
fig, ax = plt.subplots(figsize=(8, 4))
sns.heatmap(
    growth_phenotypes_to_plot.replace(value_to_int),
    cmap=cmap,
    linewidths=4,
    linecolor="white",
    annot=annot_df,
    fmt="",
    annot_kws={"fontsize": 10},
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
ax.set_yticks([i + 0.5 for i in range(len(growth_phenotypes_to_plot))])
ax.set_yticklabels(growth_phenotypes_to_plot.index, rotation=0)

# Make sure that the y-axis labels are not cut off
plt.tight_layout()

# Save the figure
plt.savefig(FILE_DIR / "pro_exomet_growth_phenotypes.png")

# Define the single color (e.g., the C-CoMP theme's green)
single_color = "#38A85E"
# Create a light palette colormap from that color
# This will create a colormap that goes from white to the specified color
green_cmap = sns.light_palette(single_color, as_cmap=True)

# Plot a heatmap of the CUE values
fig, ax = plt.subplots(figsize=(3, 4))
sns.heatmap(
    growth_phenotypes["cue"].to_frame(),
    cmap=green_cmap,
    linewidths=4,
    linecolor="white",
    fmt=".2f",
    ax=ax,
)
# Make sure that every y-tick is shown
ax.set_yticks([i + 0.5 for i in range(len(growth_phenotypes))])
ax.set_yticklabels(growth_phenotypes.index, rotation=0)
# Make sure that the y-axis labels are not cut off
plt.tight_layout()
# Save the figure
plt.savefig(FILE_DIR / "pro_exomet_cue.png")

# Plot a heatmap of the bge values
fig, ax = plt.subplots(figsize=(3, 4))
sns.heatmap(
    growth_phenotypes["bge"].to_frame(),
    cmap=green_cmap,
    linewidths=4,
    linecolor="white",
    fmt=".2f",
    ax=ax,
)
# Make sure that every y-tick is shown
ax.set_yticks([i + 0.5 for i in range(len(growth_phenotypes))])
ax.set_yticklabels(growth_phenotypes.index, rotation=0)
# Make sure that the y-axis labels are not cut off
plt.tight_layout()
# Save the figure
plt.savefig(FILE_DIR / "pro_exomet_bge.png")

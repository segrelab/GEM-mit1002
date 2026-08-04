import os
import pickle
import unittest
import warnings

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set path to the `test_files` directory
TESTFILE_DIR = os.path.join("test", "test_files")
DATA_DIR = "data"

# Define the base media
mbm_media = {
    "EX_o2_e": 20,  # O2
    "EX_na1_e": 1000,  # Na+
    "EX_cl_e": 1000,  # Cl-
    "EX_so4_e": 1000,  # Sulfate
    "EX_k_e": 1000,  # K+
    "EX_f_e": 1000,  # F-
    "EX_h2co3_e": 1000,  # H2CO3
    "EX_mg2_e": 1000,  # Mg
    "EX_ca2_e": 1000,  # Ca2+
    "EX_fe3_e": 1000,  # Fe+3
    "EX_nh3_e": 1000,  # NH3
    "EX_pi_e": 1000,  # Phosphate
    "EX_h2o_e": 1000,  # H2O
    "EX_btn_e": 1000,  # BIOT
    "EX_fol_e": 1000,  # Folate
    "EX_pydxn_e": 1000,  # yridoxol
    "EX_ribflv_e": 1000,  # Riboflavin
    "EX_thm_e": 1000,  # Thiamin
    "EX_ncam_e": 1000,  # Nicotinamide
    "EX_pnto__R_e": 1000,  # PAN
    "EX_b12_e": 1000,  # Vitamin B12
    "EX_4abz_e": 1000,  # ABEE
    "EX_cu2_e": 1000,  # Cu2+
    "EX_mn2_e": 1000,  # Mn2+
    "EX_zn2_e": 1000,  # Zn2+
    "EX_cobalt2_e": 1000,  # Co2+
}

l1_media = {
    "EX_o2_e": 20,  # O2
    "EX_na1_e": 1000,  # Na+
    "EX_no3_e": 1000,  # Nitrate
    "EX_pi_e": 1000,  # Phosphate
    "EX_fe3_e": 1000,  # Fe+3
    "EX_cl_e": 1000,  # Cl-
    "EX_mn2_e": 1000,  # Mn2+
    "EX_zn2_e": 1000,  # Zn2+
    "EX_so4_e": 1000,  # Sulfate
    "EX_cobalt2_e": 1000,  # Co2+
    "EX_cu2_e": 1000,  # Cu2+
    "EX_mobd_e": 1000,  # Molybdate
    "EX_slnt_e": 1000,  # Selenite
    "EX_ni2_e": 1000,  # Ni2+
    "EX_k_e": 1000,  # K+
    "EX_cro4_e": 1000,  # chromate
    "EX_thm_e": 1000,  # Thiamin
    "EX_btn_e": 1000,  # BIOT
    "EX_b12_e": 1000,  # Vitamin B12
    "EX_ca2_e": 1000,  # Ca2+
    "EX_mg2_e": 1000,  # Mg
}

# Load the TSV of the growth phenotypes
growth_phenotypes = pd.read_csv(
    os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"), sep="\t"
)

# Load the model
model = cobra.io.read_sbml_model("carveme_model.xml")


# Helper function for setting the media regardless if the exchange reaction is
# present in the model
# TODO: Move this to a helper file
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


# Loop through the growth phenotpes, and add the carbon source to the
# minimal media, run FBA and check if the model grows
pred_growth = []
for index, row in growth_phenotypes.iterrows():
    minimal_media = eval(row["minimal_media"]).copy()
    # Check if the model has an exchange reaction for the metabolite
    if "EX_" + str(row["bigg_id"]) + "_e" in [r.id for r in model.reactions]:
        # If it does, add the exchange reaction to the minimal media used
        minimal_media["EX_" + row["bigg_id"] + "_e"] = 1000.0
    else:
        # If it doesn't have the exchange reaction, add "No Exchange"
        pred_growth.append("No Exchange")
        # Give a warning if growth was expected
        if row["growth"] == "Yes":
            warnings.warn(
                "Model did not have an exchange reaction for "
                + row["c_source"]
                + ", but growth was expected."
            )
        continue
    # Set the media
    model.medium = clean_media(model, minimal_media)
    # Run the model
    sol = model.optimize()
    # Check if the model grows
    if sol.objective_value > 0:
        # If it does, add 'Y' to the list
        pred_growth.append("Yes")
        # Give a warning if growth was not expected
        if row["growth"] == "No":
            warnings.warn(
                "Model grew on " + row["c_source"] + ", but growth was not expected."
            )
    else:
        # If it doesn't, add 'N' to the list
        pred_growth.append("No")
        # Give a warning if growth was expected
        if row["growth"] == "Yes":
            warnings.warn(
                "Model did not grow on "
                + row["c_source"]
                + ", but growth was expected."
            )

# Add the list as a new column in the dataframe
growth_phenotypes["pred_growth"] = pred_growth

# Save the dataframe as a TSV
growth_phenotypes.to_csv(
    os.path.join(TESTFILE_DIR, "known_growth_phenotypes_w_pred.tsv"),
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

# Replace all of the "No Exchange" values with "No"
growth_phenotypes = growth_phenotypes.replace("No Exchange", "No")

# Make a dictionary for the phenotypes to numbers
value_to_int = {"Unsure": 0, "No": 1, "Yes": 2}
n = len(value_to_int)

# Make a colormap of specified colors (in numerical order for the phenotypes)
# cmap = ['gray', '#F18F01', '#399E5A'] # Gray, orange, green
cmap = ["#5E5E5E", "#FF7D0A", "#024064"]  # C-CoMP gray, orange, and dark blue

# Plot the heatmap
ax = sns.heatmap(
    growth_phenotypes.replace(value_to_int),
    cmap=cmap,
    linewidths=4,
    linecolor="white",
)

# modify colorbar:
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
    top=False,
    labeltop=True,
)

# Make sure that the y-axis labels are not cut off
plt.tight_layout()

# Save the figure
plt.savefig("exp_vs_carveme_pred_growth_phenotypes.png")

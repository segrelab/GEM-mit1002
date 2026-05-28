"""Make a scatter plot of the growth rate/yield/CUE for each single substrate
vs the abundance of that substrate in the cocktail."""

from pathlib import Path

import cobra
from gem2cue import utils
import matplotlib.pyplot as plt
import numpy as np
from numpy.f2py import main
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)


def main():
    # Load the original metabolite list to get their concentrations to sort by
    TOP_10_PATH = TOP_10_DIR / "data" / "top10_exometabolites.csv"
    top_10_exometabolites = pd.read_csv(TOP_10_PATH)
    # Filter the data to only include the Prochlorococcus marinus metabolites
    top_10_exometabolites = top_10_exometabolites[
        top_10_exometabolites["organism"] == "Prochlorococcus marinus"
    ].copy()
    # Remove the metabolites that don't support growth/that I didn't test
    # In other scripts I use the model to check this, but I don't want to have to
    # load the model here just for that, and I know from before that it's just Met
    # and Phe that don't have exchanges
    top_10_exometabolites = top_10_exometabolites[
        ~top_10_exometabolites["metabolite"].isin(["methionine", "phenylalanine"])
    ].copy()
    # Sort that list by the carbon concentration
    top_10_exometabolites = top_10_exometabolites.sort_values(
        by="carbon_concentration", ascending=False
    ).copy()

    # Load the results from the simulations from the single + cocktail simulations
    # Read the "fluxes" column as a dictionary
    singles = pd.read_csv(
        TOP_10_DIR
        / "single_and_cocktail_sims"
        / "results"
        / "single_and_cocktail_results.csv",
        converters={"fluxes": eval},
    )

    # Load the model
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    # Get the carbon exchange reactions
    c_ex_rxns = utils.get_c_ex_rxns(model)

    # Apply helper function to calcuclate CUE/GGE/BGE foreach row and merge the
    # results back into the singles DataFrame
    carbon_fates_df = singles.apply(
        lambda row: calculate_carbon_fates(row, c_ex_rxns), axis=1
    )
    singles = pd.concat([singles, carbon_fates_df], axis=1)

    # Merge the FBA results with their carbon concentrations
    singles = singles[singles["condition"] == "single"].copy()
    singles = singles.merge(
        top_10_exometabolites[["metabolite", "carbon_concentration"]],
        left_on="substrate",
        right_on="metabolite",
    )

    # Drop the "metabolite" column since it's redundant now
    singles = singles.drop(columns=["metabolite"])
    # Drop the "fluxes" column, since it is already saved elsewhere and we
    # don't need it for the plots
    singles = singles.drop(columns=["fluxes"])

    # Save the merged datframe as a csv
    singles.to_csv(
        OUT_PATH / "yields_and_carbon_concentration_per_single.csv", index=False
    )

    # Make a plot for the growth rate
    plot_y_vs_abundance(
        data=singles,
        col_name="growth_rate",
        y_label="Growth Rate on Single Substrate (1/hr)",
        title="Growth Rate vs Carbon Concentration in Cocktail",
        out_file_name="growth_rate_vs_carbon_concentration.png",
    )

    # Make a plot for the CUE
    plot_y_vs_abundance(
        data=singles,
        col_name="cue",
        y_label="Carbon Use Efficiency (CUE)",
        title="CUE vs Carbon Concentration in Cocktail",
        out_file_name="cue_vs_carbon_concentration.png",
    )

    # Make a plot for the GGE
    plot_y_vs_abundance(
        data=singles,
        col_name="gge",
        y_label="Growth Growth Efficiency (GGE)",
        title="GGE vs Carbon Concentration in Cocktail",
        out_file_name="gge_vs_carbon_concentration.png",
    )

    # Make a plot for the BGE
    plot_y_vs_abundance(
        data=singles,
        col_name="bge",
        y_label="Bacterial Growth Efficiency (BGE)",
        title="BGE vs Carbon Concentration in Cocktail",
        out_file_name="bge_vs_carbon_concentration.png",
    )


# Helper function to apply to each row
# Define a function to apply to each row
def calculate_carbon_fates(row, c_ex_rxns):
    fluxes = row["fluxes"]

    # Extract the carbon fates for the solution (both normalized and not normalized)
    c_fates = extract_c_fates_from_fluxes(
        fluxes, c_ex_rxns, co2_ex_rxn="EX_cpd00011_e0", norm=False
    )
    uptake, co2, organic_c, biomass = c_fates

    c_fates_norm = extract_c_fates_from_fluxes(
        fluxes, c_ex_rxns, co2_ex_rxn="EX_cpd00011_e0", norm=True
    )
    co2_norm, organic_c_norm, biomass_norm = c_fates_norm

    # Calculate CUE, GGE, BGE
    cue = 1 - co2 / uptake if uptake != 0 else 0
    gge = 1 - (co2 + organic_c) / uptake if uptake != 0 else 0
    bge = biomass / (biomass + co2) if (biomass + co2) != 0 else 0

    return pd.Series(
        {
            "uptake": uptake,
            "co2": co2,
            "organic_c": organic_c,
            "biomass": biomass,
            "co2_norm": co2_norm,
            "organic_c_norm": organic_c_norm,
            "biomass_norm": biomass_norm,
            "cue": cue,
            "gge": gge,
            "bge": bge,
        }
    )


def extract_c_fates_from_fluxes(fluxes, c_ex_rxns, co2_ex_rxn="EX_co2_e", norm=True):
    """This is a copy for gem2cue.utils.extract_c_fates_from_solution, but
    instead of taking a solution object, it takes a dictionary of fluxes. I
    need this because in my simulations I only save the fluxes, not the
    full solution object."""
    # Get the exchange fluxes for the current cycle
    c_ex_fluxes = {r: fluxes[r] * c for r, c in c_ex_rxns.items()}
    # Use the exchange fluxes to calculate uptake, resp, and exudation
    uptake = abs(
        sum(
            [
                flux
                for rxn, flux in c_ex_fluxes.items()
                if flux < 0 and rxn != co2_ex_rxn
            ]
        )
    )  # Should I count the co2_ex_rxn here?
    if c_ex_fluxes[co2_ex_rxn] < 0:
        # If the co2 flux is negative than the model is taking up CO2???
        co2_ex = 0
    else:
        co2_ex = c_ex_fluxes[co2_ex_rxn]
    exudation = abs(
        sum(
            [
                flux
                for rxn, flux in c_ex_fluxes.items()
                if flux > 0 and rxn != co2_ex_rxn
            ]
        )
    )
    # Calculate the biomass as everything that is not uptake or co2 release
    biomass = abs(uptake) - co2_ex - exudation
    # Normalize everything to the uptake or not
    if norm == True:
        if uptake == 0:
            co2_release_norm = 0
            exudation_norm = 0
            biomass_norm = 0
        else:
            co2_release_norm = co2_ex / uptake
            exudation_norm = exudation / uptake
            biomass_norm = biomass / uptake
        return [co2_release_norm, exudation_norm, biomass_norm]
    else:
        return [uptake, co2_ex, exudation, biomass]


def plot_y_vs_abundance(data, col_name, y_label, title, out_file_name):
    """Make a scatter plot of the growth rate/yield/CUE for each single substrate
    vs the abundance of that substrate in the cocktail."""
    # Calculate the Spearman correlation
    rho, pval = spearmanr(data["carbon_concentration"], data[col_name])

    # Make the figure
    fig, ax = plt.subplots(figsize=(8, 6))
    # Draw the scatterplot with different colors for each point/metabolite
    sns.scatterplot(
        data=data,
        x="carbon_concentration",
        y=col_name,
        hue="substrate",
        palette="tab10",
        s=100,
        edgecolor="black",
        linewidth=0.5,
    )
    # Draw a single regression line for all data (disable scatter to avaoid overlapping the colored points)
    sns.regplot(
        data=data,
        x="carbon_concentration",
        y=col_name,
        scatter=False,  # Disable the scatter points for the regression line
        line_kws={"color": "gray", "linestyle": "--"},
        ci=None,  # Don't plot the confidence interval for the regression line
    )
    # Annotate the plot with the Spearman correlation coefficient and p-value
    plt.text(
        0.05,
        0.95,
        f"Spearman ρ = {rho:.2f} (p = {pval:.2f})",
        transform=ax.transAxes,
        color="gray",
    )
    # Style
    plt.xlabel("Carbon Concentration in Cocktail (mM)", color="gray")
    plt.ylabel(y_label, color="gray")
    plt.title(title, color="gray")
    # Save
    plt.savefig(OUT_PATH / out_file_name, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()

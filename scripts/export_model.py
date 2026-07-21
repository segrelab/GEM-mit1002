import json

import cobra
import pandas as pd

# Load the model from the SBML file
model = cobra.io.read_sbml_model("model.xml")

# Export the model to JSON
cobra.io.save_json_model(model, "model.json")

# Convert to excel file
# Load the json model
with open("model.json") as f:
    model_json = json.load(f)

# Make pandas data frames for the three main components of the model
met_df = pd.DataFrame(model_json["metabolites"])
rxn_df = pd.DataFrame(model_json["reactions"])
gene_df = pd.DataFrame(model_json["genes"])


# Write a function to go from the stoichiometry dictionary to a string of the reaction with human-friendly metabolite names
# TODO: Use COBRApy's built in function for this
def build_reaction_string(stoichiometry_dict):
    reactants_side = " + ".join(
        [
            f"{abs(coeff)} {met_df.loc[met_df['id'] == met_id, 'name'].values[0]}"
            for met_id, coeff in stoichiometry_dict.items()
            if coeff < 0
        ]  # TODO: Change this into another function, so I can reuse it below
    )
    products_side = " + ".join(
        [
            f"{abs(coeff)} {met_df.loc[met_df['id'] == met_id, 'name'].values[0]}"
            for met_id, coeff in stoichiometry_dict.items()
            if coeff > 0
        ]
    )
    return f"{reactants_side} -> {products_side}"


# Apply the function to the reactions data frame and add a column for the reaction string
rxn_df["reaction"] = rxn_df["metabolites"].apply(build_reaction_string)
# Place the column third, after id and name, then keep the rest of the order the same
rxn_df = rxn_df[
    ["id", "name", "reaction"]
    + [col for col in rxn_df.columns if col not in ["id", "name", "reaction"]]
]


# Write a function to expand the annotation dictionary into sub-columns
def expand_annoation_column(df, column_name="annotation"):
    expanded = pd.json_normalize(df[column_name])
    expanded.columns = pd.MultiIndex.from_arrays(
        [[column_name] * len(expanded.columns), expanded.columns]
    )
    expanded_df = pd.concat([df.drop(columns=[column_name]), expanded], axis=1)
    return expanded_df


# Apply the function to the metabolites and reactions data frames
met_df = expand_annoation_column(met_df)
rxn_df = expand_annoation_column(rxn_df)


# Save to excel
with pd.ExcelWriter("model.xlsx") as writer:
    met_df.to_excel(writer, sheet_name="metabolites", index=False)
    rxn_df.to_excel(writer, sheet_name="reactions", index=False)
    gene_df.to_excel(writer, sheet_name="genes", index=False)

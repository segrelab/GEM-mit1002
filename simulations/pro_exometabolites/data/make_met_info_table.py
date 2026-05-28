import cobra
import pandas as pd
from pathlib import Path

FILE_DIR = Path(__file__).parent
REPO_DIR = FILE_DIR.parent.parent.parent

# Load all of the the required files
# Load the filtered data
filtered_data = pd.read_csv(FILE_DIR / "ProDiel_filtered_meanByTimepoint.csv")
# Load my mapping of met IDs to names
met_mapping = pd.read_csv(FILE_DIR.parent / "metabolite_id_map.csv")
# Load the model
model = cobra.io.read_sbml_model(REPO_DIR / "model.xml")

# Make a new df, with the column "name" being all of the unique names in the filtered data
met_info_table = pd.DataFrame()
met_info_table["name"] = filtered_data["CleanName"].drop_duplicates()

# For each name, get the corresponding met ID from the mapping file, and add it to the table
met_info_table["id"] = met_info_table["name"].map(
    lambda x: met_mapping.loc[met_mapping["name"] == x, "id"].values[0]
)

# For each met ID, check if the "EX_" + met ID reaction is in the model, and add a column "EX_in_model" with True or False
met_info_table["EX_in_model"] = met_info_table["id"].apply(
    lambda x: True if f"EX_{x}_e0" in [r.id for r in model.reactions] else False
)

# For each met ID, check if the c0 version of the met ID is in the model, and add a column "c0_in_model" with True or False
met_info_table["c0_in_model"] = met_info_table["id"].apply(
    lambda x: True if f"{x}_c0" in [m.id for m in model.metabolites] else False
)


# For each met ID, add the metabolite's formula from the model, if it exists, and add it to a column "formula"
def get_formula(met_id):
    c0_id = f"{met_id}_c0"
    if c0_id in [m.id for m in model.metabolites]:
        return model.metabolites.get_by_id(c0_id).formula
    else:
        return None


met_info_table["formula"] = met_info_table["id"].apply(get_formula)

# Print the markdown table
print(met_info_table.to_markdown(index=False))

# Save the table as a csv
met_info_table.to_csv(FILE_DIR / "met_info_table.csv", index=False)

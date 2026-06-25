import sys
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Import plot_styles.py from the root of the repo
sys.path.append(str(REPO_ROOT))
from plot_styles import summer_colors

# Load the results
df = pd.read_csv(OUT_PATH / "results.csv", header=[0, 1], index_col=0)

# Select the desired columns (e.g. CUE, BGE, growth rate)
cue_df = df.xs("cue", axis=1, level=1)  # Get just the CUE columns

# Remove any columns with all NaN values
cue_df = cue_df.dropna(axis=1, how="all")

# Melt the datafrom into long format for categorical plotting
cue_df_melted = cue_df.reset_index().melt(id_vars="index", value_name="CUE")

# TODO: Sort the biomass components by chemical class? By stoichiometric coefficient?

# Define a custom heatmap
# Define the colors
low_color = summer_colors["pink"]  # Color for the low extreme
mid_color = summer_colors["dark_tan"]  # Color for the midpoint
high_color = summer_colors["teal"]  # Color for the high extreme

# Define the relative positions of the colors (0 = min, 0.5 = midpoint, 1 = max)
colors = [(0, low_color), (0.5, mid_color), (1, high_color)]

# Create the custom colormap
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", colors)

# Create the plot
g = sns.catplot(
    data=cue_df_melted,
    x="component",
    y="CUE",
    hue="index",
    palette=custom_cmap,
    height=6,
    aspect=2,
    legend_out=True,
)
# TODO: Make the x-tick labels use the human-friendly names for the metabolites?
# Turn the x-tick labels
g.set_xticklabels(rotation=90)
g.ax.set_xlabel("Biomass Component", color="gray")
g.ax.set_ylabel("CUE (mmol C/ mmol C)", color="gray")  # TODO: Check the units
# set_plot_style(g.ax)

# Shrink the plot to make room for the x labels
g.fig.subplots_adjust(bottom=0.3)

# TODO: Give the legend a better title- show the linear values, make it continuous?

# Save the plot
plt.savefig(OUT_PATH / "cue_vs_biomass_perturbation.png")

import os
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))

# Add the project root to the system path to import my plot_styles file
sys.path.append(PROJECT_ROOT)
# Import the plot style from the plot_styles.py file
from plot_styles import summer_colors, set_plot_style

# Load the data
data = pd.read_csv(os.path.join(FILE_DIR, "growth_match_summary.csv"))

# Skip rows that say "Error" in the % Match column
data = data[data["% Match"] != "ERROR"]

# Convert columns to numeric type
# TODO: Do I need any other columns numerical? Are there any that can't be?
cols_to_fix = ["Matches", "% Match", "PR Number", "Unbounded Flux Reactions"]
data[cols_to_fix] = data[cols_to_fix].apply(pd.to_numeric, errors='coerce')

# Sort the data by PR number and reset the index
data = data.sort_values("PR Number").reset_index(drop=True)

# Pick one color per series so the line and its axis can be matched by eye
matches_color = summer_colors["teal"]
flux_color = summer_colors["pink"]

# Create a figure with twin y axes
fig, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()
# Scale the unbounded flux axis to show detail near zero
ax2.set_yscale("symlog", linthresh=1)

# Plot the number of matches on the left axis
ax1.plot(
    data.index,
    data["Matches"],
    marker="o",
    markersize=4,
    linestyle="-",
    color=matches_color,
)

# Plot the number of arbitrarily large reactions on the right axis
ax2.plot(
    data.index,
    data["Unbounded Flux Reactions"],
    marker="o",
    markersize=4,
    linestyle="-",
    color=flux_color,
)

# Make it so it looks like the matches line is "on top"
# Move ax1 to a higher z-order than ax2
ax1.set_zorder(ax2.get_zorder() + 1)
# Make ax1's background transparent so ax2 is still visible behind it
ax1.patch.set_visible(False)

# Apply the shared style (gray bottom axis, no top/right spines, gray text)
set_plot_style(ax1)

# Titles and labels (set on ax1 so set_plot_style's gray text applies)
ax1.set_title("Model Performance Over Time")
ax1.set_xlabel("Pull Request Number")
ax1.set_ylabel("Growth Phenotypes Matching Experimental Data")
ax1.set_ylim(0, 50)  # See the full range
ax2.set_ylabel("Unique Reactions with Flux > 100 (Log Scale)")

# Thin out the x-tick labels: with ~120 points, labeling every PR is unreadable,
# so show every Nth PR number instead
step = max(1, len(data) // 15)
tick_positions = data.index[::step]
ax1.set_xticks(tick_positions)
ax1.set_xticklabels(data["PR Number"].iloc[::step], rotation=45, ha="right")

# Color the left axis to match the matches line
ax1.spines["left"].set_color(matches_color)
ax1.tick_params(axis="y", colors=matches_color)
ax1.yaxis.label.set_color(matches_color)

# Color the right axis to match the unbounded flux line.
# set_plot_style hides ax2's spines, and its right spine is the one we want, so
# style ax2 by hand instead of calling set_plot_style on it.
ax2.spines["top"].set_visible(False)
ax2.spines["left"].set_visible(False)
ax2.spines["right"].set_color(flux_color)
ax2.tick_params(axis="y", colors=flux_color)
ax2.yaxis.label.set_color(flux_color)

# Save the plot
fig.tight_layout()
fig.savefig(
    os.path.join(FILE_DIR, "match_over_time.png"),
    dpi=150,
)

import os
import sys
import warnings

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
from adjustText import adjust_text

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))

# Add the project root to the system path to import my plot_styles file
sys.path.append(PROJECT_ROOT)
# Import the plot style from the plot_styles.py file
from plot_styles import set_plot_style, summer_colors

# Define which PRs to highlight on the plot
# These are PRs that caused a significant change in the number of matches or unbounded flux reactions
GROWTH_HIGHLIGHT_PRS = {
    160: "Coneect dead-ends to glyoxylate",
    192: "Removed duplicate thiamine reaction",
    281: "Standardize amino acid transporters (no longer grows on leucine and isoleucine)",
    309: "Add Na+ symporters for leucine and isoleucine",
}
FLUX_HIGHLIGHT_PRS = {
    211: "Remove reactions with glucose anomers",
    344: "Fixed the GAPDHN infeasible loop",
}

# Load the data
data = pd.read_csv(os.path.join(FILE_DIR, "growth_match_summary.csv"))

# Skip rows that say "Error" in the % Match column
data = data[data["% Match"] != "ERROR"]

# Convert columns to numeric type
# TODO: Do I need any other columns numerical? Are there any that can't be?
cols_to_fix = ["Matches", "% Match", "PR Number", "Unbounded Flux Reactions"]
data[cols_to_fix] = data[cols_to_fix].apply(pd.to_numeric, errors="coerce")

# Sort the data by PR number and reset the index
data = data.sort_values("PR Number").reset_index(drop=True)

# Pick one color per series so the line and its axis can be matched by eye
matches_color = summer_colors["teal"]
flux_color = summer_colors["pink"]

# Create a figure with twin y axes
fig, ax1 = plt.subplots(figsize=(8, 5))
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
ax1.set_ylim(0, 55)  # See the full range
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

# ---- Annotate the highlighted PRs -------------------------------------------
# The x-axis is the row index, not the PR number, so map PR number -> x position
pr_to_index = dict(zip(data["PR Number"], data.index))


def darken(color, factor=0.7):
    """Return a darker shade of a color so the star stands out from its line."""
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)


# Give the top of the growth axis a little headroom for labels
ax1.set_ylim(0, 58)
# Finalize the layout first so the data->display transforms used below (to put
# the flux labels into ax1's coordinate system) are correct.
fig.tight_layout()
fig.canvas.draw()

# Collect all highlight labels in ONE coordinate system (ax1) so adjustText can
# place growth and flux labels together without them colliding. Stars are drawn
# on each series' own axis; the flux point position is transformed into ax1
# coords only so its label and leader line land in the right place.
label_texts = []


def add_highlights(source_ax, y_col, highlights, series_color):
    star_color = darken(series_color)
    for pr in highlights:
        if pr not in pr_to_index:
            warnings.warn(f"Highlight PR #{pr} is not in the plotted data; skipping.")
            continue
        x = pr_to_index[pr]
        y = data.loc[x, y_col]
        # Bigger star, darker shade of the line color, white halo to separate it
        source_ax.plot(
            x,
            y,
            marker="*",
            markersize=22,
            color=star_color,
            markeredgecolor="white",
            markeredgewidth=1.4,
            zorder=7,
        )
        # Position the label in ax1 data coords (transform if it's a flux point)
        disp = source_ax.transData.transform((x, y))
        x1, y1 = ax1.transData.inverted().transform(disp)
        label_texts.append(
            ax1.text(
                x1,
                y1,
                f"PR {pr}",
                color=star_color,
                fontsize=12,
                fontweight="bold",
                ha="center",
                va="center",
                zorder=8,
            )
        )


add_highlights(ax1, "Matches", GROWTH_HIGHLIGHT_PRS, matches_color)
add_highlights(ax2, "Unbounded Flux Reactions", FLUX_HIGHLIGHT_PRS, flux_color)

# Points for the labels to avoid: every plotted marker of both series, all in
# ax1 coords so labels don't sit on either line.
avoid_x, avoid_y = [], []
for col, src in (("Matches", ax1), ("Unbounded Flux Reactions", ax2)):
    for x, y in zip(data.index, data[col]):
        dx, dy = ax1.transData.inverted().transform(src.transData.transform((x, y)))
        avoid_x.append(dx)
        avoid_y.append(dy)

adjust_text(
    label_texts,
    x=avoid_x,
    y=avoid_y,
    ax=ax1,
    arrowprops=dict(arrowstyle="-", color="0.5", lw=0.8),
    expand=(1.4, 1.8),
    force_text=(0.4, 0.7),
    ensure_inside_axes=True,
    min_arrow_len=6,
)

# Save the plot
fig.savefig(os.path.join(FILE_DIR, "match_over_time.png"), dpi=300)

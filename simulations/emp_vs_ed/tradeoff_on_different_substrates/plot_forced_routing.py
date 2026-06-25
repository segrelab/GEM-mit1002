import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[2]
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Import plot_styles.py from the root of the repo
sys.path.append(str(REPO_ROOT))
from plot_styles import summer_colors

# Set some parameters for the plot
width = 0.35
emp_colors = summer_colors["teal"]
ed_colors = summer_colors["pink"]

# Load the results
results_df = pd.read_csv(OUT_PATH / "forced_routing_results.csv", index_col=0)

# Extract the data fromt the datafreame as lists
substrates = results_df.index.tolist()
emp_vals = results_df["emp"].tolist()
ed_vals = results_df["ed"].tolist()

# Get the x positions for the bars
x = np.arange(len(substrates))

# Plot the data
fig, ax = plt.subplots(figsize=(10, 6))
bars_emp = ax.bar(x - width / 2, emp_vals, width, label="EMP", color=emp_colors)
bars_ed = ax.bar(x + width / 2, ed_vals, width, label="ED", color=ed_colors)

# Plot style
ax.set_xlabel("Substrate")
ax.set_ylabel("Growth rate (h$^{-1}$)")
ax.set_title("Growth rate under forced pathway routing")
ax.set_xticks(x)
ax.set_xticklabels(substrates, rotation=45, ha="right")
ax.legend()
fig.tight_layout()

# Save the figure
fig.savefig(OUT_PATH / "forced_routing_plot.png", dpi=150)

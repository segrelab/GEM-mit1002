import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy import stats

matplotlib.rcParams.update(
    {
        "font.size": 11,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]
IN_PATH = FILE_PATH / "results"
OUT_PATH = FILE_PATH / "figures"
OUT_PATH.mkdir(exist_ok=True)

# Import plot_styles.py from the root of the repo
sys.path.append(str(REPO_ROOT))
from plot_styles import summer_colors

# Anchor colours drawn from the "Summer I Turned Pretty" palette (plus a muted
# mauve), ordered to flow as a gradient. build_palette() uses these directly
# when there are few categories, and interpolates between them when there are
# more categories than anchors (e.g. one colour per substrate).
PALETTE_ANCHORS = [
    summer_colors["teal"],
    summer_colors["light_blue"],
    summer_colors["green"],
    summer_colors["yellow"],
    summer_colors["pink"],
    summer_colors["dark_pink"],
    "#9B6A9E",  # muted mauve
    summer_colors["dark_tan"],
]


def build_palette(categories) -> dict:
    """Map an ordered list of categories to colours from the summer palette.

    With no more categories than anchor colours, the anchors are used directly;
    with more, colours are interpolated along the palette so every category
    gets a distinct hue."""
    cats = list(dict.fromkeys(categories))  # unique, order-preserving
    n = len(cats)
    if n == 0:
        return {}
    if n <= len(PALETTE_ANCHORS):
        cols = PALETTE_ANCHORS[:n]
    else:
        cmap = LinearSegmentedColormap.from_list("summer", PALETTE_ANCHORS)
        cols = [cmap(i / (n - 1)) for i in range(n)]
    return dict(zip(cats, cols))


# Fixed order for the entry-point categories, with colours from the palette.
ENTRY_POINT_ORDER = [
    "Glycolysis",
    "Pyruvate",
    "Acetyl-CoA",
    "TCA — α-KG",
    "TCA — OAA",
    "TCA — Succinyl-CoA",
    "Aromatic catabolism",
]
ENTRY_POINT_COLORS = build_palette(ENTRY_POINT_ORDER)


def main():
    # Load the substrate panel
    substrate_df = pd.read_csv(IN_PATH / "substrate_panel.csv")
    # Load the Growth and CUE results
    summary_df = pd.read_csv(IN_PATH / "growth_and_cue.csv", index_col=0)
    # Rename the "name" column in the substrate_df to "substrate" so that it can be merged with the results_df
    substrate_df = substrate_df.rename(columns={"name": "substrate"})
    # Merge the results with the substrate info
    merged_df = pd.merge(summary_df, substrate_df, on="met_id")

    # Extract only the data for the anchor level
    anchor_level = sorted(merged_df["o2_level"].unique(), reverse=True)[0]
    anchor_df = merged_df[merged_df["o2_level"] == anchor_level]

    # Plot the bar charts of growth rate and CUE, coloured by entry point into central metabolism
    # Only for the "anchor level" (unlimited O2)
    plot_growth_cue(anchor_df, OUT_PATH, "growth_and_cue")

    # Plot the correlation between growth rate and CUE for the anchor level
    plot_growth_cue_correlation(
        anchor_df, OUT_PATH, "growth_vs_cue_scatter_anchor", metric="cue"
    )
    # Plot the correlation between growth rate and BGE for the anchor level
    plot_growth_cue_correlation(
        anchor_df, OUT_PATH, "growth_vs_bge_scatter_anchor", metric="bge"
    )

    # Plot the correlations across all oxygen levels
    # CUE, colored by entry point
    plot_growth_cue_correlation(
        merged_df, OUT_PATH, "growth_vs_cue_scatter", metric="cue"
    )
    # CUE, coloured by substrate instead of entry point
    plot_growth_cue_correlation(
        merged_df,
        OUT_PATH,
        "growth_vs_cue_scatter_by_substrate",
        metric="cue",
        color_by="substrate",
    )
    # BGE, colored by entry point
    plot_growth_cue_correlation(
        merged_df, OUT_PATH, "growth_vs_bge_scatter", metric="bge"
    )
    # BGE, coloured by substrate instead of entry point
    plot_growth_cue_correlation(
        merged_df,
        OUT_PATH,
        "growth_vs_bge_scatter_by_substrate",
        metric="bge",
        color_by="substrate",
    )


def plot_growth_cue(summary_df: pd.DataFrame, out_path: Path, filename: str) -> None:
    df = summary_df.sort_values("growth_rate", ascending=False).reset_index()
    colors = [ENTRY_POINT_COLORS.get(c, "#9E9E9E") for c in df["entry_point"]]
    x = np.arange(len(df))

    fig, (ax_g, ax_c) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

    ax_g.bar(x, df["growth_rate"], color=colors, edgecolor="white", linewidth=0.6)
    ax_g.set_ylabel("Growth rate (h⁻¹)")
    ax_g.set_title(
        "MIT1002 growth rate and carbon-use efficiency across substrates",
        fontsize=12,
        pad=8,
    )
    ax_g.margins(x=0.01)

    ax_c.bar(x, df["cue"], color=colors, edgecolor="white", linewidth=0.6)
    ax_c.set_ylabel("Carbon-use efficiency")
    ax_c.set_ylim(0, 1)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(df["substrate"], rotation=45, ha="right", fontsize=9)

    handles = [
        mpatches.Patch(facecolor=ENTRY_POINT_COLORS[lbl], edgecolor="white", label=lbl)
        for lbl in ENTRY_POINT_ORDER
        if lbl in df["entry_point"].values
    ]
    ax_g.legend(
        handles=handles,
        title="Entry point into\ncentral metabolism",
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
        title_fontsize=9.5,
    )

    fig.tight_layout()
    fig.savefig(out_path / f"{filename}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_growth_cue_correlation(
    summary_df: pd.DataFrame,
    out_path: Path,
    filename: str,
    metric="cue",
    color_by="entry_point",
) -> None:
    df = summary_df.reset_index()

    # Check that the metric column exists
    if metric not in df.columns:
        raise ValueError(f"Column '{metric}' not found in the DataFrame")

    # Check that the color_by column exists
    if color_by not in df.columns:
        raise ValueError(f"Column '{color_by}' not found in the DataFrame")

    # Define a label for the y-axis based on the metric
    if metric == "cue":
        y_label = "Carbon Use efficiency"
    elif metric == "bge":
        y_label = "Bacterial Growth Efficiency"
    else:
        y_label = metric

    # Oxygen levels, highest first. The dot is anchored at the highest level
    # (where growth and CUE are ~perfectly correlated, i.e. along the diagonal)
    # and arrows trace how each substrate moves as O2 drops to the next level.
    levels = sorted(df["o2_level"].unique(), reverse=True)
    anchor_level = levels[0]

    anchor = df[df["o2_level"] == anchor_level]

    # Build the colour mapping for whichever column we're colouring by.
    # entry_point uses a fixed category order; substrate is ordered by anchor
    # growth rate so the legend reads top-to-bottom like the diagonal.
    if color_by == "entry_point":
        categories = [c for c in ENTRY_POINT_ORDER if c in df[color_by].values]
        legend_title = "Entry point into\ncentral metabolism"
    else:
        categories = anchor.sort_values("growth_rate", ascending=False)[
            color_by
        ].tolist()
        legend_title = color_by.replace("_", " ").capitalize()
    palette = build_palette(categories)

    def cat_color(value):
        return palette.get(value, "#9E9E9E")

    colors = [cat_color(c) for c in anchor[color_by]]
    x = anchor["growth_rate"].values.astype(float)
    y = anchor[metric].values.astype(float)

    # Correlation reported at the anchor (highest-O2) level.
    r, p = stats.pearsonr(x, y)
    rho, p_s = stats.spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    # Diagonal dots: one per substrate at the highest O2 level.
    ax.scatter(x, y, c=colors, s=120, edgecolor="white", linewidth=0.6, zorder=4)

    # Least-squares fit line through the anchor-level points.
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, slope * xs + intercept, color="#555555", ls="--", lw=1.2, zorder=2)

    # For each substrate, draw arrows from one O2 level to the next one down,
    # following the trajectory: anchor_level -> ... -> lowest level.
    for substrate, sdf in df.groupby("substrate"):
        sdf = sdf.set_index("o2_level")
        color = cat_color(sdf[color_by].iloc[0])
        for lo_hi, lo_lo in zip(levels[:-1], levels[1:]):
            if lo_hi not in sdf.index or lo_lo not in sdf.index:
                continue
            x0, y0 = sdf.loc[lo_hi, "growth_rate"], sdf.loc[lo_hi, metric]
            x1, y1 = sdf.loc[lo_lo, "growth_rate"], sdf.loc[lo_lo, metric]
            ax.annotate(
                "",
                xy=(x1, y1),
                xytext=(x0, y0),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=color,
                    lw=1.3,
                    alpha=0.8,
                    shrinkA=4,
                    shrinkB=2,
                ),
                zorder=3,
            )

    # Substrate labels at the anchor (diagonal) dots. Skipped when colouring by
    # substrate, since the legend already identifies each point. Default offset
    # is to the right; near-coincident points get a manual override.
    if color_by != "substrate":
        label_offsets = {  # substrate: (dx, dy, ha, va)
            "Galactose": (0, 9, "center", "bottom"),  # sits almost on Glucose
        }
        for _, row in anchor.iterrows():
            dx, dy, ha, va = label_offsets.get(
                row["substrate"], (5, 0, "left", "center")
            )
            ax.annotate(
                row["substrate"],
                (row["growth_rate"], row[metric]),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=7,
                ha=ha,
                va=va,
                color="#555555",
            )

    # Annotations don't trigger autoscaling, so expand the axes to cover all
    # O2 levels (the arrow endpoints), with a small margin.
    gx = df["growth_rate"].astype(float)
    gy = df[metric].astype(float)
    xpad = 0.05 * (gx.max() - gx.min())
    ypad = 0.05 * (gy.max() - gy.min())
    ax.set_xlim(gx.min() - xpad, gx.max() + xpad)
    ax.set_ylim(gy.min() - ypad, gy.max() + ypad)

    ax.set_xlabel("Growth rate (h⁻¹)")
    ax.set_ylabel(y_label)
    ax.set_title(f"Growth rate vs. {y_label} across substrates", fontsize=12, pad=8)

    # Correlation metrics (at the highest O2 level).
    txt = (
        f"At O₂ = {anchor_level:g}\n"
        f"Pearson  r = {r:.2f}  (p = {p:.2g})\n"
        f"Spearman ρ = {rho:.2f}  (p = {p_s:.2g})"
    )
    # Placed inside the bottom-right corner (clear of the up-left arrow fan) so
    # it doesn't collide with the legend, which can be long when colouring by
    # substrate.
    ax.text(
        0.98,
        0.03,
        txt,
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="right",
        linespacing=1.4,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc", alpha=0.9),
    )

    # Legend: one entry per category, plus the fit line and the dot/arrow keys.
    handles = [
        mpatches.Patch(facecolor=cat_color(lbl), edgecolor="white", label=lbl)
        for lbl in categories
    ]
    handles.append(
        Line2D([0], [0], color="#555555", ls="--", lw=1.2, label="linear fit")
    )
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="#777777",
            lw=0,
            markersize=7,
            label=f"dot: O₂ = {anchor_level:g}",
        )
    )
    handles.append(
        Line2D(
            [0],
            [0],
            color="#777777",
            lw=1.3,
            label="arrow: decreasing O₂",
        )
    )
    # Shrink the legend text when there are many categories (e.g. per substrate).
    leg_fontsize = 9 if len(categories) <= 10 else 7
    ax.legend(
        handles=handles,
        title=legend_title,
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=leg_fontsize,
        title_fontsize=9.5,
    )

    fig.tight_layout()
    # Save the figure as a PNG
    fig.savefig(out_path / f"{filename}.png", dpi=300, bbox_inches="tight")
    # Save the figure as an SVG
    fig.savefig(out_path / f"{filename}.svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}  (Pearson r={r:.2f}, p={p:.2g})")


if __name__ == "__main__":
    main()

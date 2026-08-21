import ast
import re
import sys
from pathlib import Path

import cobra
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
FATES_OUT_PATH = OUT_PATH / "carbon_fates"
FATES_OUT_PATH.mkdir(exist_ok=True)
MODEL_PATH = REPO_ROOT / "model.xml"

# Organic byproducts whose max carbon flux (mmol C/gDW/h) across the *whole*
# dataset is below this are lumped into a single grey "Other" segment in the
# carbon fates bars, keeping trace byproducts out of the legend/palette.
BYPRODUCT_FLUX_THRESHOLD = 1.0

# Import the shared plot styles from tools/
sys.path.append(str(REPO_ROOT))
from tools.plot_styles import (
    carbon_fates_bar,
    ccomp_colors,
    set_plot_style,
    summer_colors,
)

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


def build_palette(categories, anchors=PALETTE_ANCHORS) -> dict:
    """Map an ordered list of categories to colours from a palette.

    With no more categories than anchor colours, the anchors are used directly;
    with more, colours are interpolated along the palette so every category
    gets a distinct hue."""
    cats = list(dict.fromkeys(categories))  # unique, order-preserving
    n = len(cats)
    if n == 0:
        return {}
    if n <= len(anchors):
        cols = anchors[:n]
    else:
        cmap = LinearSegmentedColormap.from_list("palette", anchors)
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

# Anchors for the organic-byproduct palette, deliberately excluding
# summer_colors["teal"] and summer_colors["yellow"] since carbon_fates_bar
# hard-codes those for the "Biomass" and "CO2" segments.
BYPRODUCT_PALETTE_ANCHORS = [
    summer_colors["light_blue"],
    summer_colors["green"],
    summer_colors["pink"],
    summer_colors["dark_pink"],
    "#9B6A9E",  # muted mauve
    summer_colors["dark_tan"],
    ccomp_colors["orange"],
    ccomp_colors["dark_blue"],
]
OTHER_BYPRODUCT_COLOR = "#B0B0B0"


def build_ex_name_map(model) -> dict:
    """Map exchange reaction id -> readable metabolite name.

    Mirrors the renaming done in run_growth_cue.build_exchange_df, so
    byproducts are labelled the same way here as in exchange_fluxes.csv."""
    names = {}
    for r in model.reactions:
        if not r.id.startswith("EX_"):
            continue
        met_name = next(iter(r.metabolites)).name
        suffix = " [e0]"
        if met_name.endswith(suffix):
            met_name = met_name[: -len(suffix)]
        names[r.id] = met_name
    return names


def parse_c_ex_fluxes(raw) -> dict:
    """Parse the "c_ex_fluxes" column back into a {rxn_id: flux} dict.

    It's written by run_growth_cue.py as the repr() of a plain {str: float}
    dict. Older result files wrapped the values in np.float64(...), which
    isn't valid for ast.literal_eval, so strip that wrapper if present."""
    if pd.isna(raw):
        return {}
    cleaned = re.sub(r"np\.float64\(([^)]+)\)", r"\1", raw)
    return ast.literal_eval(cleaned)


def byproduct_flux_df(df: pd.DataFrame, ex_name_map: dict) -> pd.DataFrame:
    """Expand the "c_ex_fluxes" column into one column per byproduct name.

    Values are carbon flux (mmol C/gDW/h); rows/index match `df`."""
    rows = []
    for raw in df["c_ex_fluxes"]:
        named = {}
        for rxn_id, flux in parse_c_ex_fluxes(raw).items():
            name = ex_name_map.get(rxn_id, rxn_id)
            named[name] = named.get(name, 0.0) + flux
        rows.append(named)
    return pd.DataFrame(rows, index=df.index).fillna(0.0)


def build_byproduct_palette(
    summary_df: pd.DataFrame, ex_name_map: dict, threshold=BYPRODUCT_FLUX_THRESHOLD
) -> dict:
    """Build a {byproduct_name: color} palette covering the whole dataset.

    Byproducts are ranked by their max carbon flux across every substrate and
    O2 level so the most prominent ones get distinct colors; used everywhere
    a carbon fates plot is made so a given byproduct (e.g. Acetate) always
    gets the same color, no matter which substrate's plot it appears in.
    Byproducts that never exceed `threshold` share a single grey "Other"
    color instead of eating into the distinct-color budget."""
    all_byproducts = byproduct_flux_df(summary_df, ex_name_map)
    max_flux = all_byproducts.abs().max().sort_values(ascending=False)
    major = max_flux[max_flux >= threshold].index.tolist()
    palette = build_palette(major, anchors=BYPRODUCT_PALETTE_ANCHORS)
    palette["Other"] = OTHER_BYPRODUCT_COLOR
    return palette


def carbon_fates_byproduct_columns(
    sub_df: pd.DataFrame, ex_name_map: dict, palette: dict
) -> pd.DataFrame:
    """Per-byproduct carbon flux columns for one substrate's carbon fates bar.

    Byproducts with their own entry in `palette` keep their own column;
    anything else is lumped into "Other". All-zero columns are dropped so
    the legend only lists byproducts this substrate actually releases."""
    byp_df = byproduct_flux_df(sub_df, ex_name_map)
    # Order by rank in the (already flux-ranked) palette, not by first-seen
    # order in the data, so segment order matches across every substrate
    named_cols = [c for c in palette if c != "Other" and c in byp_df.columns]
    other_cols = [c for c in byp_df.columns if c not in palette]
    out = byp_df[named_cols].copy()
    if other_cols:
        out["Other"] = byp_df[other_cols].sum(axis=1)
    return out.loc[:, out.abs().sum(axis=0) > 1e-9]


def main():
    # Load the substrate panel
    substrate_df = pd.read_csv(IN_PATH / "substrate_panel.csv")
    # Load the Growth and CUE results
    summary_df = pd.read_csv(IN_PATH / "growth_and_cue.csv", index_col=0)
    # Rename the "name" column in the substrate_df to "substrate" so that it can be merged with the results_df
    substrate_df = substrate_df.rename(columns={"name": "substrate"})
    # Merge the results with the substrate info
    merged_df = pd.merge(summary_df, substrate_df, on="met_id")
    # Round the O2 percent values to 2 decimal places
    merged_df["o2_percent"] = merged_df["o2_percent"].round(2)

    # Load the model (needed to turn exchange reaction ids in "c_ex_fluxes"
    # into readable metabolite names) and build a palette for the organic
    # byproducts, shared across every carbon fates plot so a given byproduct
    # (e.g. Acetate) is always the same color
    print("Loading model for byproduct names...")
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    ex_name_map = build_ex_name_map(model)
    byproduct_palette = build_byproduct_palette(merged_df, ex_name_map)

    # Extract only the data for the anchor level
    anchor_level = sorted(merged_df["o2_bound"].unique(), reverse=True)[0]
    anchor_df = merged_df[merged_df["o2_bound"] == anchor_level]

    # Extract the data for the "pre-set" oxygen levels
    # Assume that the "pre-set" oxygen levels are those where the bound is a round number (e.g. 10, 20, 50)
    # So its modulo 10 is 0
    # And add the o2_bound == 1 level
    pre_set_df = merged_df[merged_df["o2_bound"].mod(10) == 0]
    pre_set_df = pd.concat([pre_set_df, merged_df[merged_df["o2_bound"] == 1]])

    # Extract the data for the "percentile" oxygen levels
    # Assume that the "percentile" oxygen levels are those where the percent is a round number (e.g. 10, 20, 30, 40, 50)
    # So its modulo 10 is 0
    # Round the percent value first in case things are very slightly off
    percentile_df = merged_df[round(merged_df["o2_percent"], 2).mod(10) == 0]

    # Plot the bar charts of growth rate and CUE, coloured by entry point into central metabolism
    # Only for the "anchor level" (unlimited O2)
    plot_growth_cue(anchor_df, OUT_PATH, "growth_and_cue")

    # Make carbon fate plots for each substrate
    for met in pre_set_df["substrate"].unique():
        # Subset the data for the current substrate
        met_df = pre_set_df[pre_set_df["substrate"] == met]
        # Set the index to be the o2_bound column
        met_df = met_df.set_index("o2_bound")
        # Expand the organic byproducts into one column per compound
        byp_df = carbon_fates_byproduct_columns(met_df, ex_name_map, byproduct_palette)
        # Rename the carbon fates columns to be what the plotting function expects
        met_df.rename(
            columns={
                "biomass_c": "biomass",
                "co2_flux": "co2",
            },
            inplace=True,
        )
        # Subset the columns needed for the plot
        met_df = pd.concat([met_df[["biomass"]], byp_df, met_df[["co2"]]], axis=1)
        # Make the plot
        g = carbon_fates_bar(met_df, byproduct_colors=byproduct_palette)
        # Set the title to be something specific to the metabolite
        g.set_title(f"Carbon Fates for {met}")
        # Set the x axis label to be more descriptive
        g.set_xlabel("Oxygen Uptake Bound")
        # Set the y axis label
        g.set_ylabel("Carbon Flux (mmol C/gDW/h)")
        # Save the plot
        plt.tight_layout()
        plt.savefig(
            FATES_OUT_PATH / f"carbon_fates_{met}_set_o2_bounds.png",
            dpi=300,
            bbox_inches="tight",
        )
    # Do the same for the percentile O2 data
    for met in percentile_df["substrate"].unique():
        # Subset the data for the current substrate
        met_df = percentile_df[percentile_df["substrate"] == met]
        # Set the index to be the o2_percent column
        met_df = met_df.set_index("o2_percent")
        # Order to be in descenting o2_percent
        met_df = met_df.sort_index(ascending=False)
        # Expand the organic byproducts into one column per compound
        byp_df = carbon_fates_byproduct_columns(met_df, ex_name_map, byproduct_palette)
        # Rename the carbon fates columns to be what the plotting function expects
        met_df.rename(
            columns={
                "biomass_c": "biomass",
                "co2_flux": "co2",
            },
            inplace=True,
        )
        # Subset the columns needed for the plot
        met_df = pd.concat([met_df[["biomass"]], byp_df, met_df[["co2"]]], axis=1)
        # Make the plot
        g = carbon_fates_bar(met_df, byproduct_colors=byproduct_palette)
        # Set the title to be something specific to the metabolite
        g.set_title(f"Carbon Fates for {met}")
        # Set the x axis label to be more descriptive
        g.set_xlabel("Percent Oxygen Saturation")
        # Set the y axis label
        g.set_ylabel("Carbon Flux (mmol C/gDW/h)")
        # Save the plot
        plt.tight_layout()
        plt.savefig(
            FATES_OUT_PATH / f"carbon_fates_{met}_o2_percentiles.png",
            dpi=300,
            bbox_inches="tight",
        )

    # Plot the correlation between growth rate and CUE for the anchor level
    plot_growth_cue_correlation(
        anchor_df, OUT_PATH, "growth_vs_cue_scatter_anchor", metric="cue"
    )
    # Plot the correlation between growth rate and BGE for the anchor level
    plot_growth_cue_correlation(
        anchor_df, OUT_PATH, "growth_vs_bge_scatter_anchor", metric="bge"
    )
    # Plot the correlation between growth rate and GGE for the anchor level
    plot_growth_cue_correlation(
        anchor_df, OUT_PATH, "growth_vs_gge_scatter_anchor", metric="gge"
    )

    # Plot the correlations across all pre-set oxygen levels
    # CUE, colored by entry point
    plot_growth_cue_correlation(
        pre_set_df, OUT_PATH, "growth_vs_cue_scatter", metric="cue"
    )
    # CUE, coloured by substrate instead of entry point
    plot_growth_cue_correlation(
        pre_set_df,
        OUT_PATH,
        "growth_vs_cue_scatter_by_substrate",
        metric="cue",
        color_by="substrate",
    )
    # BGE, colored by entry point
    plot_growth_cue_correlation(
        pre_set_df, OUT_PATH, "growth_vs_bge_scatter", metric="bge"
    )
    # BGE, coloured by substrate instead of entry point
    plot_growth_cue_correlation(
        pre_set_df,
        OUT_PATH,
        "growth_vs_bge_scatter_by_substrate",
        metric="bge",
        color_by="substrate",
    )
    # GGE, colored by entry point
    plot_growth_cue_correlation(
        pre_set_df, OUT_PATH, "growth_vs_gge_scatter", metric="gge"
    )
    # GGE, coloured by substrate instead of entry point
    plot_growth_cue_correlation(
        pre_set_df,
        OUT_PATH,
        "growth_vs_gge_scatter_by_substrate",
        metric="gge",
        color_by="substrate",
    )

    # Plot the correlations across all percentile oxygen levels
    # CUE, coloured by substrate
    plot_growth_cue_correlation(
        percentile_df,
        OUT_PATH,
        "growth_vs_cue_scatter_by_substrate_percentiles",
        metric="cue",
        o2_level_col="o2_percent",
        color_by="substrate",
    )
    # BGE, coloured by substrate
    plot_growth_cue_correlation(
        percentile_df,
        OUT_PATH,
        "growth_vs_bge_scatter_by_substrate_percentiles",
        metric="bge",
        o2_level_col="o2_percent",
        color_by="substrate",
    )
    # GGE, coloured by substrate
    plot_growth_cue_correlation(
        percentile_df,
        OUT_PATH,
        "growth_vs_gge_scatter_by_substrate_percentiles",
        metric="gge",
        o2_level_col="o2_percent",
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
    o2_level_col="o2_bound",
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
        y_label = "Carbon Use Efficiency"
    elif metric == "bge":
        y_label = "Bacterial Growth Efficiency"
    elif metric == "gge":
        y_label = "Gross Growth Efficiency"
    else:
        y_label = metric

    # Define a label for O2 levels used
    if o2_level_col == "o2_bound":
        o2_label = "O2 Lower Bound"
    elif o2_level_col == "o2_percent":
        o2_label = "O2 Saturation Percentile"
    else:
        o2_label = o2_level_col

    # Find the "anchor" level- the highest O2 level to plot the dots on the diagonal
    # If you set O2 bounds, all will have data for the same "o2_bound" value
    # But if you are doing percentiles, use the o2_percent column
    levels = sorted(df[o2_level_col].unique(), reverse=True)
    anchor_level = levels[0]
    for met in df["met_id"].unique():
        if anchor_level not in df[df["met_id"] == met][o2_level_col].values:
            raise ValueError(
                f"Metabolite '{met}' is missing data for the anchor level '{anchor_level}'"
            )
    # Extract the data for the anchor level
    anchor = df[df[o2_level_col] == anchor_level]

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
        sdf = sdf.set_index(o2_level_col)
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
        f"At {o2_label} = {anchor_level:g}\n"
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
            label=f"dot: {o2_label} = {anchor_level:g}",
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

    # Set the plot style (gray axes, ticks, labels, title and legend)
    set_plot_style(ax)

    # Save the figure as a PNG
    fig.savefig(out_path / f"{filename}.png", dpi=300, bbox_inches="tight")
    # Save the figure as an SVG
    fig.savefig(out_path / f"{filename}.svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {filename}.png/svg (Pearson r={r:.2f}, p={p:.2g})")


if __name__ == "__main__":
    main()

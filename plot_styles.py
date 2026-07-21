import matplotlib.pyplot as plt

# Define the colors from the C-CoMP pallette
ccomp_colors = {
    "light_blue": "#3CB3C0",
    "dark_blue": "#024064",
    "orange": "#FF6C2C",
    "light_orange": "#FFBB62",  # Not technically a C-CoMP color, but I used it on my ISME poster
}

# Define the colors for the "Summer I Turned Pretty" palette
summer_colors = {
    "teal": "#5B8C8F",
    "light_blue": "#BBD5E9",
    "green": "#7E9B6B",
    "pink": "#E57E66",
    "dark_pink": "#893625",
    "yellow": "#EBB309",
    "dark_tan": "#C9BC9B",
    "light_tan": "#FBF9EA",
}


# Define the style for the plots (gray axes, no top or right axis lines)
def set_plot_style(g):
    # Make the axis lines gray
    g.spines["bottom"].set_color("gray")
    g.spines["left"].set_color("gray")
    # Make the tick marks gray
    g.tick_params(axis="x", colors="gray")
    g.tick_params(axis="y", colors="gray")
    # Remove the top and right axis lines
    g.spines["top"].set_visible(False)
    g.spines["right"].set_visible(False)
    # Make all text (axis labels, tick labels, title, and legend) gray
    g.xaxis.label.set_color("gray")
    g.yaxis.label.set_color("gray")
    g.title.set_color("gray")
    if g.get_legend() is not None:
        for text in g.get_legend().get_texts():
            text.set_color("gray")


def carbon_fates_bar(data, byproduct_colors=None):
    """Stacked bar plot of the fate of carbon: biomass, organic byproducts, CO2.

    By default `data` must have exactly the columns "biomass", "organic_c",
    and "co2", and the organic byproducts are shown as a single lumped bar.

    Pass `byproduct_colors` (an ordered {column_name: color} dict, e.g. from
    `build_byproduct_palette`) to instead split that bar into one segment per
    organic byproduct. `data` must then have a "biomass" and "co2" column
    plus one column per byproduct named in `byproduct_colors` (any columns
    named in `byproduct_colors` but absent from `data` are simply skipped, so
    the same palette can be reused across substrates that release different
    byproducts). Passing the same `byproduct_colors` dict across multiple
    calls keeps a given byproduct's color consistent across figures.
    """
    if byproduct_colors is not None:
        # Check that the required columns are present
        assert "biomass" in data.columns and "co2" in data.columns
        byproduct_cols = [c for c in byproduct_colors if c in data.columns]
        # Check that there aren't any unexpected byproduct columns
        extra_cols = set(data.columns) - set(byproduct_cols) - {"biomass", "co2"}
        assert not extra_cols, (
            f"Columns not found in byproduct_colors palette: {extra_cols}"
        )
        # Set the column order
        data = data[["biomass"] + byproduct_cols + ["co2"]]
        colors = (
            [summer_colors["teal"]]
            + [byproduct_colors[c] for c in byproduct_cols]
            + [summer_colors["yellow"]]
        )
        custom_labels = ["Biomass"] + byproduct_cols + ["CO2"]
    else:
        # Check that the column names are correct
        assert set(data.columns) == set(["co2", "organic_c", "biomass"])
        # Set the column order
        data = data[["biomass", "organic_c", "co2"]]
        colors = [
            summer_colors["teal"],
            summer_colors["light_blue"],
            summer_colors["yellow"],
        ]
        custom_labels = ["Biomass", "Organic C", "CO2"]
    # Plot the stacked bar plot
    g = data.plot(
        kind="bar",
        stacked=True,
        color=colors,
    )
    # Move the legend outside of the plot
    lgd = plt.legend(
        bbox_to_anchor=(1.25, 0.5),
        loc="center right",
        borderaxespad=0.0,
        ncol=1,
        labels=custom_labels,
    )
    # Make the legend text gray too
    for text in lgd.get_texts():
        text.set_color("gray")
    # Adjust the bottom margin
    plt.subplots_adjust(bottom=0.2)
    # Style
    set_plot_style(g)
    # Title the plot and make it gray
    g.set_title("Fate of Carbon", color="gray")

    # Return the plot
    return g

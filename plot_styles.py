import matplotlib.pyplot as plt

# Define the colors from the C-CoMP pallette
ccomp_colors = {
    "light_blue": "#3CB3C0",
    "dark_blue": "#024064",
    "orange": "#FF6C2C",
    "light_orange": "#FFBB62",  # Not technically a C-CoMP color, but I used it on my ISME poster
}

# Define the colors for the "Summer I Turned Pretty" palette
tsitp_colors = {
    "dark_blue": "#2D5F80",
    "med_blue": "#4D82A6",
    "light_blue": "#BBD5E9",
    "dark_green": "#5A8C72",
    "light_green": "#BDDACC",
    "dark_orange": "#C0604A",
    "light_orange": "#E57E66",
    "dark_pink": "#C4788A",
    "light_pink": "#F2BFB5",
    "dark_yellow": "#EBB309",
    "light_yellow": "#F6D5A6",
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


def carbon_fates_bar(data):
    # Check that the column names are correct
    assert set(data.columns) == set(["co2", "organic_c", "biomass"])
    # Set the column order
    data = data[["biomass", "organic_c", "co2"]]
    # Plot the stacked bar plot
    g = data.plot(
        kind="bar",
        stacked=True,
        color=[
            ccomp_colors["dark_blue"],
            ccomp_colors["light_blue"],
            ccomp_colors["light_orange"],
        ],
    )
    # Move the legend outside of the plot
    custom_labels = ["Biomass", "Organic C", "CO2"]
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

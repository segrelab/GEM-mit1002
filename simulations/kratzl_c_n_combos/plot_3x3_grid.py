"""Nitrogen source x carbon source grid: experiment vs. FBA prediction.

The Kratzl nitrogen-source screen is the only genuine factorial block in
``data/known_growth_phenotypes.tsv``: three nitrogen sources crossed with three
carbon sources, all in ``marine_broth_wo_yeast_and_peptone_no_n``. That is the
one place a 2D grid earns its keep, because conditions failing for a single
shared reason line up as a whole row or column instead of scattering through a
61-row list.

Every cell is split on the diagonal. The upper-left triangle is the
experiment, the lower-right triangle is the model, so a cell reading as one
solid block is an agreement and a two-tone cell is a mismatch. Mismatched
cells also carry a heavy outline: colour is already carrying the data, so the
annotation gets its own channel.

Cells with no exchange reaction in the model are hatched rather than drawn as
a no-growth prediction. The distinction matters for what this figure is
arguing. A metabolite the model cannot import produces zero growth for a
trivial reason, and scoring that as a correct "No" would inflate specificity
and let a reconstruction gap masquerade as a limit of constraint-based
modelling.

Deliberately absent: any mapping of growth *rate* to colour. Against binary
data the magnitude carries no validation signal and invites reading a
precision that is not there. Predictions within 10x of the growth threshold
get a small dot instead, since those are the calls most likely to flip on an
unrelated curation change.

Run from the repository root::

    python simulations/kratzl_c_n_combos/plot_3x3_grid.py
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch, Polygon, Rectangle

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
sys.path.insert(0, PROJECT_ROOT)

from tools.phenotypes import (  # noqa: E402
    C_SOURCE_IDS,
    GROWTH_THRESHOLD,
    N_SOURCE_IDS,
    evaluate_phenotypes,
    format_summary,
    load_phenotypes,
    summarise,
)
from tools.plot_styles import summer_colors  # noqa: E402

FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
MODEL_PATH = os.path.join(PROJECT_ROOT, "model.xml")

#: The medium this screen was run in. Both the carbon and the nitrogen source
#: come from the condition rather than the base medium.
N_SCREEN_MEDIUM = "marine_broth_wo_yeast_and_peptone_no_n"

#: Set to a letter (e.g. "A") when this grid is assembled into a multi-panel
#: figure. Left off for the standalone version, where a lone panel letter with
#: no siblings just raises questions.
PANEL_LABEL = None

GROWTH_FILL = summer_colors["teal"]
#: Solid rather than the paler light_tan, so that "did not grow" and "hatched
#: because the model has no exchange reaction" separate by fill as well as by
#: texture. Three states need three clearly different values.
NO_GROWTH_FILL = summer_colors["dark_tan"]
CELL_EDGE = summer_colors["dark_tan"]
HATCH_EDGE = summer_colors["dark_tan"]
DISCORDANT_EDGE = summer_colors["dark_pink"]
TEXT_COLOR = "#333333"

#: Grid cell size, in inches. Only nine cells carry the result, so they get to
#: be large; the key and legend are explanatory furniture and stay subordinate.
CELL = 0.70
PAD = 0.08
#: Width of the white stripe along the diagonal, in cell units. Only needs to
#: be wide enough to keep the split legible when both halves agree; at this
#: cell size anything larger reads as a design element rather than a divider.
GAP = 0.03

#: Size of the key's example cell as a fraction of a grid cell. Kept well under
#: 1.0 so the example reads as an explanation rather than as a tenth data
#: point competing with the nine real ones.
KEY_CELL_SCALE = 0.55

#: The key's outline, slightly lighter than the grid's so it does not turn into
#: a thick frame at the smaller cell size.
KEY_DISCORDANT_LW = 1.6

#: Mark predictions that land within NEAR_THRESHOLD_FACTOR of the growth
#: threshold. Only ever drawn on cells whose prediction is actually shown, so
#: it stays off the hatched "no exchange reaction" cells. Set False to drop the
#: markers entirely.
SHOW_NEAR_THRESHOLD = True

#: Room to reserve for the key's two annotation labels, in inches.
KEY_TEXT_LEFT = 0.68
KEY_TEXT_RIGHT = 0.42

#: Vertical space between the annotated example cell and the colour key. They
#: are one explanatory unit, so this is deliberately tight.
KEY_LEGEND_GAP = 0.12

LABEL_SIZE = 8
HEADER_SIZE = 8
PANEL_SIZE = 9
DISCORDANT_LW = 2.0


def _axes_inches(fig, x, y, width, height):
    """Add an axes positioned in inches from the top-left of the figure.

    Sizing the axes explicitly rather than letting the layout engine do it is
    what keeps the cells square and identically sized.
    """
    fig_width, fig_height = fig.get_size_inches()
    return fig.add_axes(
        [
            x / fig_width,
            1 - (y + height) / fig_height,
            width / fig_width,
            height / fig_height,
        ]
    )


def _half_style(state, representable=True):
    """Patch keywords for one triangle of a cell."""
    if not representable or state is None:
        return dict(
            facecolor="white", hatch="////", edgecolor=HATCH_EDGE, linewidth=0.5
        )
    if state == "Yes":
        return dict(facecolor=GROWTH_FILL, edgecolor="none", linewidth=0.0)
    if state == "No":
        return dict(facecolor=NO_GROWTH_FILL, edgecolor="none", linewidth=0.0)
    return dict(facecolor="white", hatch="....", edgecolor=HATCH_EDGE, linewidth=0.5)


def _draw_cell(ax, col, row, record):
    """Draw one split cell. ``record`` is a row of the evaluated table."""
    x, y = col, row
    representable = record["predicted"] is not None

    upper_left = [
        (x + PAD, y + PAD),
        (x + 1 - PAD - GAP, y + PAD),
        (x + PAD, y + 1 - PAD - GAP),
    ]
    lower_right = [
        (x + 1 - PAD, y + PAD + GAP),
        (x + 1 - PAD, y + 1 - PAD),
        (x + PAD + GAP, y + 1 - PAD),
    ]

    ax.add_patch(
        Polygon(upper_left, closed=True, zorder=2, **_half_style(record["growth"]))
    )
    ax.add_patch(
        Polygon(
            lower_right,
            closed=True,
            zorder=2,
            **_half_style(record["predicted"], representable),
        )
    )

    if record["discordant"]:
        ax.add_patch(
            Rectangle(
                (x + PAD - 0.03, y + PAD - 0.03),
                1 - 2 * PAD + 0.06,
                1 - 2 * PAD + 0.06,
                fill=False,
                edgecolor=DISCORDANT_EDGE,
                linewidth=DISCORDANT_LW,
                zorder=4,
            )
        )

    # Skipped when the cell is hatched: flagging a prediction as marginal makes
    # no sense on a condition the model cannot represent, where any nonzero rate
    # comes from the parts of the medium that did get added.
    if record["near_threshold"] and representable and SHOW_NEAR_THRESHOLD:
        # Outlined rather than filled, so it stays visible on either fill.
        ax.plot(
            [x + 0.72],
            [y + 0.72],
            marker="o",
            markersize=2.6,
            markerfacecolor="white",
            markeredgecolor=CELL_EDGE,
            markeredgewidth=0.5,
            zorder=5,
        )


def _draw_matrix(ax, cells, row_labels, col_labels):
    """Draw the matrix of split cells. ``cells`` maps (row, col) -> record."""
    n_rows, n_cols = len(row_labels), len(col_labels)

    for (row, col), record in cells.items():
        _draw_cell(ax, col, row, record)

    # Untested combinations stay empty but keep a hairline, so "not tested" is
    # distinguishable from "tested and did not grow".
    for row in range(n_rows):
        for col in range(n_cols):
            if (row, col) not in cells:
                ax.add_patch(
                    Rectangle(
                        (col + PAD, row + PAD),
                        1 - 2 * PAD,
                        1 - 2 * PAD,
                        fill=False,
                        edgecolor=CELL_EDGE,
                        linewidth=0.4,
                        linestyle=(0, (1, 1.6)),
                        zorder=1,
                    )
                )

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.invert_yaxis()
    ax.set_xticks([c + 0.5 for c in range(n_cols)])
    ax.set_yticks([r + 0.5 for r in range(n_rows)])
    ax.set_yticklabels(row_labels, fontsize=LABEL_SIZE, color=TEXT_COLOR)
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(
        col_labels,
        fontsize=HEADER_SIZE,
        color=TEXT_COLOR,
        rotation=40,
        ha="left",
        rotation_mode="anchor",
    )
    ax.tick_params(length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)


def build_grid(results):
    """Map the nitrogen-source screen onto a nitrogen x carbon grid."""
    rows = list(N_SOURCE_IDS)
    cols = list(C_SOURCE_IDS)

    cells = {}
    for _, record in results.iterrows():
        met_ids = set(record["met_id"])
        for row_index, n_id in enumerate(rows):
            for col_index, c_id in enumerate(cols):
                if met_ids == {n_id, c_id}:
                    cells[(row_index, col_index)] = record

    return cells, [N_SOURCE_IDS[n] for n in rows], [C_SOURCE_IDS[c] for c in cols]


def _draw_key(ax, width, height, cell):
    """A single annotated split cell explaining which half is which.

    Data coordinates here are inches, so ``cell`` sets the example cell's true
    printed size and it can be kept smaller than a grid cell.
    """
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.invert_yaxis()
    ax.axis("off")

    gap = GAP * cell
    x = KEY_TEXT_LEFT + 0.14
    y = 0.09

    ax.add_patch(
        Polygon(
            [(x, y), (x + cell - gap, y), (x, y + cell - gap)],
            closed=True,
            facecolor=GROWTH_FILL,
            zorder=2,
        )
    )
    ax.add_patch(
        Polygon(
            [(x + cell, y + gap), (x + cell, y + cell), (x + gap, y + cell)],
            closed=True,
            facecolor=NO_GROWTH_FILL,
            edgecolor="none",
            zorder=2,
        )
    )
    ax.add_patch(
        Rectangle(
            (x - 0.02, y - 0.02),
            cell + 0.04,
            cell + 0.04,
            fill=False,
            edgecolor=DISCORDANT_EDGE,
            linewidth=KEY_DISCORDANT_LW,
            zorder=3,
        )
    )
    ax.annotate(
        "experiment",
        xy=(x + 0.20 * cell, y + 0.18 * cell),
        xytext=(x - 0.11, y - 0.01),
        fontsize=LABEL_SIZE,
        color=TEXT_COLOR,
        ha="right",
        va="center",
        arrowprops=dict(arrowstyle="-", lw=0.6, color=CELL_EDGE),
    )
    ax.annotate(
        "model",
        xy=(x + 0.80 * cell, y + 0.82 * cell),
        xytext=(x + cell + 0.11, y + cell + 0.01),
        fontsize=LABEL_SIZE,
        color=TEXT_COLOR,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", lw=0.6, color=CELL_EDGE),
    )


def _legend_entries(cells):
    """Legend handles for the states that actually occur in the panel.

    Built from the data rather than hard-coded so the legend never explains a
    symbol the reader cannot find. With only the nitrogen screen plotted, that
    drops the "uncertain" and "not tested" entries.
    """
    experimental = {record["growth"] for record in cells.values()}
    predicted = {record["predicted"] for record in cells.values()}
    n_slots = len(N_SOURCE_IDS) * len(C_SOURCE_IDS)

    entries = [
        (
            "Yes" in experimental or "Yes" in predicted,
            Patch(facecolor=GROWTH_FILL, label="Growth"),
        ),
        (
            "No" in experimental or "No" in predicted,
            Patch(facecolor=NO_GROWTH_FILL, label="No growth"),
        ),
        (
            "Unsure" in experimental,
            Patch(
                facecolor="white",
                hatch="....",
                edgecolor=HATCH_EDGE,
                label="Growth uncertain (experiment)",
            ),
        ),
        (
            None in predicted,
            Patch(
                facecolor="white",
                hatch="////",
                edgecolor=HATCH_EDGE,
                label="No exchange reaction in model",
            ),
        ),
        (
            any(record["discordant"] for record in cells.values()),
            Patch(
                facecolor="none",
                edgecolor=DISCORDANT_EDGE,
                linewidth=1.7,
                label="Experiment and model disagree",
            ),
        ),
        (
            len(cells) < n_slots,
            Patch(
                facecolor="none",
                edgecolor=CELL_EDGE,
                linestyle=(0, (1, 1.6)),
                label="Not tested",
            ),
        ),
    ]
    return [handle for include, handle in entries if include]


def _draw_legend(ax, handles):
    ax.axis("off")
    ax.legend(
        handles=handles,
        loc="upper left",
        frameon=False,
        fontsize=LABEL_SIZE,
        handlelength=1.3,
        handleheight=1.3,
        labelspacing=0.55,
        borderpad=0,
    )


def generate_grid_figure(model, output_stem="kratzl_c_n_grid"):
    # Only the nine conditions in this screen are solved. There is no reason to
    # run the other 52 phenotypes to draw a 3x3 grid.
    phenotypes = load_phenotypes()
    subset = phenotypes[phenotypes["minimal_media"] == N_SCREEN_MEDIUM]
    results = evaluate_phenotypes(model, phenotypes=subset)
    summary = summarise(results)
    print(format_summary(summary))

    cells, row_labels, col_labels = build_grid(results)
    handles = _legend_entries(cells)

    margin = 0.18
    label_width = 0.80
    # Just enough to clear the rotated column labels. Over-reserving here only
    # creates whitespace that bbox_inches="tight" crops back off, which throws
    # off anything positioned relative to the figure height.
    header = 0.60
    key_cell = CELL * KEY_CELL_SCALE
    key_width = KEY_TEXT_LEFT + 0.14 + key_cell + 0.14 + KEY_TEXT_RIGHT
    key_height = key_cell + 0.20
    legend_width = 2.95
    # Matches what the legend actually occupies at this font size and
    # spacing; an over-estimate leaves the axes reserving blank space that
    # bbox_inches="tight" then keeps.
    legend_height = 0.06 + 0.22 * len(handles)

    matrix_x = margin + label_width
    matrix_width = len(col_labels) * CELL
    matrix_height = len(row_labels) * CELL
    matrix_y = header + 0.28

    # The example cell sits this far into its own axes, because the
    # "experiment" label is drawn to its left. Backing the axes off by the same
    # amount puts the example cell on the same vertical line as the legend
    # swatches, so everything in the right-hand column aligns.
    key_indent = KEY_TEXT_LEFT + 0.14
    right_x = matrix_x + matrix_width + 0.90
    key_x = right_x - key_indent

    # Centre the key and legend together on the grid rows. Measuring against
    # the rows rather than the rows plus column labels keeps this independent
    # of how tall the rotated labels happen to render.
    matrix_bottom = matrix_y + matrix_height
    block_height = key_height + KEY_LEGEND_GAP + legend_height
    key_y = matrix_y + (matrix_height - block_height) / 2
    legend_y = key_y + key_height + KEY_LEGEND_GAP

    fig_width = right_x + max(key_width - key_indent, legend_width) + margin
    fig_height = max(matrix_bottom, legend_y + legend_height) + margin + 0.10
    fig = plt.figure(figsize=(fig_width, fig_height))

    ax_matrix = _axes_inches(fig, matrix_x, matrix_y, matrix_width, matrix_height)
    _draw_matrix(ax_matrix, cells, row_labels, col_labels)

    if PANEL_LABEL:
        fig.text(
            margin / fig_width,
            1 - (matrix_y - header - 0.24) / fig_height,
            PANEL_LABEL,
            fontsize=PANEL_SIZE,
            fontweight="bold",
            va="top",
            ha="left",
            color=TEXT_COLOR,
        )

    _draw_key(
        _axes_inches(fig, key_x, key_y, key_width, key_height),
        key_width,
        key_height,
        key_cell,
    )
    _draw_legend(
        _axes_inches(fig, right_x, legend_y, legend_width, legend_height), handles
    )

    os.makedirs(FIGURES_DIR, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(
            os.path.join(FIGURES_DIR, f"{output_stem}.{extension}"),
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)

    # The classified table backs the numbers quoted in the caption.
    results.drop(columns=["met_id"]).to_csv(
        os.path.join(FIGURES_DIR, f"{output_stem}_classified.tsv"),
        sep="\t",
        index=False,
    )
    with open(os.path.join(FIGURES_DIR, f"{output_stem}_summary.txt"), "w") as handle:
        handle.write(format_summary(summary) + "\n")
        handle.write(f"Growth threshold: {GROWTH_THRESHOLD} 1/hr\n")

    return results, summary


def main():
    import cobra

    model = cobra.io.read_sbml_model(MODEL_PATH)
    results, _ = generate_grid_figure(model)

    discordant = results[results["discordant"]]
    if not discordant.empty:
        print("\nMismatches:")
        print(
            discordant[["c_source", "growth", "predicted", "category"]].to_string(
                index=False
            )
        )

    unrepresentable = results[results["category"] == "no_exchange"]
    if not unrepresentable.empty:
        print("\nNot representable (no exchange reaction):")
        print(unrepresentable[["c_source", "missing_exchanges"]].to_string(index=False))


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    main()

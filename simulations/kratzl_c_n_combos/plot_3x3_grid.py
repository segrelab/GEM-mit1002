"""Manuscript figure: experimental vs. FBA-predicted growth phenotypes.

Every cell is split on the diagonal. The upper-left triangle is the
experiment, the lower-right triangle is the model, so a cell that reads as one
solid block is an agreement and a two-tone cell is a mismatch. Mismatched
cells additionally carry a heavy outline, because colour alone is doing enough
work already and the discordant cells are the point of the figure.

Three panels, one per experimental design, all sharing one encoding:

A. Nitrogen source x carbon source (Kratzl). The only genuine factorial block
   in the dataset, and the only layout where a 2D grid earns its keep: cells
   that fail for one shared reason line up as a row or column.
B. Amino acid with and without pyruvate (Kratzl), in a medium that already
   contains ammonium and nitrate, so this is a carbon-source screen.
C. Single-substrate carbon screens, as metabolite x medium. Using medium as
   the second axis (rather than co-substrate) surfaces the cases where the
   same metabolite was scored differently by different labs.

Cells with no exchange reaction in the model are hatched rather than drawn as
a no-growth prediction. That distinction matters for the argument the figure
is making: a missing transporter is an incomplete reconstruction, not a
limitation of constraint-based modelling, and the two should not share a
colour.

Deliberately absent: any mapping of growth *rate* to colour. Against binary
data the magnitude carries no validation signal and invites reading a
precision that is not there. Marginal predictions (within 10x of the growth
threshold) get a small dot instead, since those are the calls that will flip
on an unrelated curation change.

Run from the repository root::

    python -m scripts.generate_phenotype_figure
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch, Polygon, Rectangle

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from tools.phenotypes import (  # noqa: E402
    C_SOURCE_IDS,
    GROWTH_THRESHOLD,
    N_SOURCE_IDS,
    evaluate_phenotypes,
    format_summary,
    summarise,
)
from tools.plot_styles import ccomp_colors  # noqa: E402

RESULTS_DIR = os.path.join(PROJECT_ROOT, "scripts", "results")
MODEL_PATH = os.path.join(PROJECT_ROOT, "model.xml")

PYRUVATE = "cpd00020"
N_SCREEN_MEDIUM = "marine_broth_wo_yeast_and_peptone_no_n"
C_SCREEN_MEDIUM = "marine_broth_wo_yeast_and_peptone"

#: Media shown in panel C, in display order.
SINGLE_SUBSTRATE_MEDIA = {
    "mbm": "MBM",
    "l1": "L1",
    "promm_no_c": "ProMM",
    "swm": "SWM",
}

#: Shorter labels for the handful of names that would dominate the axis.
SHORT_NAMES = {
    "DHPS (dihydroxypropanesulfonate)": "DHPS",
    "3-methyl-2-oxobutanoic acid": "3-methyl-2-oxobutanoate",
    "3-methyl-2-oxopentanoic acid": "3-methyl-2-oxopentanoate",
    "4-methyl-2-oxopentanoic acid": "4-methyl-2-oxopentanoate",
    "4-hydroxybenzoic acid": "4-hydroxybenzoate",
    "Glycerol-3-phosphate": "Glycerol-3-P",
}

GROWTH_FILL = ccomp_colors["dark_blue"]
NO_GROWTH_FILL = "#E8E6DF"
HATCH_EDGE = "#8C8C8C"
DISCORDANT_EDGE = ccomp_colors["orange"]
CELL_EDGE = "#C9C6BC"

#: Panel C has one row per metabolite, which is much taller than it is wide.
#: Wrapping it into this many side-by-side blocks keeps the figure a usable
#: shape for a journal page.
PANEL_C_BLOCKS = 2

CELL = 0.30
PAD = 0.08
GAP = 0.05
LABEL_SIZE = 8
HEADER_SIZE = 8
PANEL_SIZE = 9


def short_name(name: str) -> str:
    return SHORT_NAMES.get(name, name)


def _axes_inches(fig, x, y, width, height):
    """Add an axes positioned in inches from the top-left of the figure.

    Sizing the axes rather than letting the layout engine do it is what keeps
    every cell in every panel exactly the same size.
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
            facecolor="white", hatch="////", edgecolor=HATCH_EDGE, linewidth=0.4
        )
    if state == "Yes":
        return dict(facecolor=GROWTH_FILL, edgecolor="none", linewidth=0.0)
    if state == "No":
        return dict(facecolor=NO_GROWTH_FILL, edgecolor=CELL_EDGE, linewidth=0.4)
    return dict(facecolor="white", hatch="....", edgecolor=HATCH_EDGE, linewidth=0.4)


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
                linewidth=1.7,
                zorder=4,
            )
        )

    if record["near_threshold"]:
        ax.plot(
            [x + 0.72],
            [y + 0.72],
            marker="o",
            markersize=2.2,
            color="white",
            zorder=5,
        )


def _draw_matrix(ax, cells, row_labels, col_labels, rotate_columns=False):
    """Draw a matrix of split cells. ``cells`` maps (row, col) -> record."""
    n_rows, n_cols = len(row_labels), len(col_labels)

    for (row, col), record in cells.items():
        _draw_cell(ax, col, row, record)

    # Cells with no experiment at all stay empty, but get a hairline so the
    # reader can tell "not tested" from "tested and no growth".
    for row in range(n_rows):
        for col in range(n_cols):
            if (row, col) not in cells:
                ax.add_patch(
                    Rectangle(
                        (col + PAD, row + PAD),
                        1 - 2 * PAD,
                        1 - 2 * PAD,
                        fill=False,
                        edgecolor="#E0DED7",
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
    ax.set_yticklabels(row_labels, fontsize=LABEL_SIZE, color="#333333")
    ax.xaxis.set_ticks_position("top")
    if rotate_columns:
        ax.set_xticklabels(
            col_labels,
            fontsize=HEADER_SIZE,
            color="#333333",
            rotation=40,
            ha="left",
            rotation_mode="anchor",
        )
    else:
        ax.set_xticklabels(col_labels, fontsize=HEADER_SIZE, color="#333333")
    ax.tick_params(length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _figure_label(fig, x, y, text):
    """Panel letter, positioned in inches from the top-left of the figure.

    Placed on the figure rather than the axes so it clears the column headers,
    which sit above the axes and vary in height between panels.
    """
    fig_width, fig_height = fig.get_size_inches()
    fig.text(
        x / fig_width,
        1 - y / fig_height,
        text,
        fontsize=PANEL_SIZE,
        fontweight="bold",
        va="top",
        ha="left",
        color="#333333",
    )


def _split_blocks(cells, row_labels, n_blocks):
    """Break a tall matrix into side-by-side blocks of rows.

    Panel C has one row per metabolite, which is far taller than it is wide.
    Wrapping it into blocks keeps the figure a usable shape without changing
    what is plotted.
    """
    if n_blocks < 2:
        return [(cells, row_labels)]

    chunk = -(-len(row_labels) // n_blocks)
    blocks = []
    for index in range(n_blocks):
        offset = index * chunk
        rows = row_labels[offset : offset + chunk]
        if not rows:
            continue
        blocks.append(
            (
                {
                    (row - offset, col): record
                    for (row, col), record in cells.items()
                    if offset <= row < offset + len(rows)
                },
                rows,
            )
        )
    return blocks


def build_nitrogen_carbon_panel(results):
    """Panel A: nitrogen source x carbon source."""
    rows = list(N_SOURCE_IDS)
    cols = list(C_SOURCE_IDS)
    subset = results[results["minimal_media"] == N_SCREEN_MEDIUM]

    cells = {}
    for _, record in subset.iterrows():
        met_ids = set(record["met_id"])
        for row_index, n_id in enumerate(rows):
            for col_index, c_id in enumerate(cols):
                if met_ids == {n_id, c_id}:
                    cells[(row_index, col_index)] = record
    return cells, [N_SOURCE_IDS[n] for n in rows], [C_SOURCE_IDS[c] for c in cols]


def build_pyruvate_panel(results):
    """Panel B: amino acid alone vs. amino acid plus pyruvate."""
    subset = results[results["minimal_media"] == C_SCREEN_MEDIUM]

    records = {}
    for _, record in subset.iterrows():
        met_ids = list(record["met_id"])
        base = [m for m in met_ids if m != PYRUVATE]
        if len(base) != 1:
            continue
        name = short_name(record["c_source"].split(",")[0].strip())
        column = 1 if PYRUVATE in met_ids else 0
        records[(name, column)] = record

    row_labels = sorted({name for name, _ in records})
    cells = {
        (row_labels.index(name), column): record
        for (name, column), record in records.items()
    }
    return cells, row_labels, ["alone", "+ pyruvate"]


def build_medium_panel(results):
    """Panel C: single-substrate carbon screens, metabolite x medium."""
    subset = results[
        results["minimal_media"].isin(SINGLE_SUBSTRATE_MEDIA)
        & (results["met_id"].apply(len) == 1)
    ]

    records = {}
    for _, record in subset.iterrows():
        name = short_name(record["c_source"])
        records[(name, record["minimal_media"])] = record

    # Group by experimental outcome, then alphabetical within group. Sorting
    # discordant cells to the top would make the figure look arranged.
    def group(name):
        calls = {
            record["growth"]
            for (metabolite, _), record in records.items()
            if metabolite == name
        }
        if "Yes" in calls:
            return 0
        if "No" in calls:
            return 1
        return 2

    names = sorted({name for name, _ in records})
    row_labels = sorted(names, key=lambda name: (group(name), name.lower()))
    col_keys = list(SINGLE_SUBSTRATE_MEDIA)

    cells = {
        (row_labels.index(metabolite), col_keys.index(medium)): record
        for (metabolite, medium), record in records.items()
    }
    return cells, row_labels, [SINGLE_SUBSTRATE_MEDIA[k] for k in col_keys]


def _draw_key(ax):
    """A single annotated split cell explaining which half is which."""
    ax.set_xlim(0, 3.4)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()
    ax.axis("off")

    x, y = 0.9, 0.05
    size = 0.9
    upper_left = [
        (x, y),
        (x + size - GAP, y),
        (x, y + size - GAP),
    ]
    lower_right = [
        (x + size, y + GAP),
        (x + size, y + size),
        (x + GAP, y + size),
    ]
    ax.add_patch(Polygon(upper_left, closed=True, facecolor=GROWTH_FILL, zorder=2))
    ax.add_patch(
        Polygon(
            lower_right,
            closed=True,
            facecolor=NO_GROWTH_FILL,
            edgecolor=CELL_EDGE,
            linewidth=0.4,
            zorder=2,
        )
    )
    ax.add_patch(
        Rectangle(
            (x - 0.03, y - 0.03),
            size + 0.06,
            size + 0.06,
            fill=False,
            edgecolor=DISCORDANT_EDGE,
            linewidth=1.7,
            zorder=3,
        )
    )
    ax.annotate(
        "experiment",
        xy=(x + 0.18, y + 0.16),
        xytext=(x - 0.25, y - 0.02),
        fontsize=LABEL_SIZE,
        color="#333333",
        ha="right",
        va="center",
        arrowprops=dict(arrowstyle="-", lw=0.6, color="#8C8C8C"),
    )
    ax.annotate(
        "model",
        xy=(x + size - 0.18, y + size - 0.16),
        xytext=(x + size + 0.25, y + size + 0.02),
        fontsize=LABEL_SIZE,
        color="#333333",
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", lw=0.6, color="#8C8C8C"),
    )


def _draw_legend(ax):
    ax.axis("off")
    handles = [
        Patch(facecolor=GROWTH_FILL, label="Growth"),
        Patch(facecolor=NO_GROWTH_FILL, edgecolor=CELL_EDGE, label="No growth"),
        Patch(
            facecolor="white",
            hatch="....",
            edgecolor=HATCH_EDGE,
            label="Growth uncertain (experiment)",
        ),
        Patch(
            facecolor="white",
            hatch="////",
            edgecolor=HATCH_EDGE,
            label="No exchange reaction in model",
        ),
        Patch(
            facecolor="none",
            edgecolor=DISCORDANT_EDGE,
            linewidth=1.7,
            label="Experiment and model disagree",
        ),
        Patch(
            facecolor="none",
            edgecolor="#E0DED7",
            linestyle=(0, (1, 1.6)),
            label="Not tested",
        ),
    ]
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


def generate_phenotype_figure(model, output_stem="exp_vs_pred_phenotype_matrix"):
    results = evaluate_phenotypes(model)
    summary = summarise(results)
    print(format_summary(summary))

    a_cells, a_rows, a_cols = build_nitrogen_carbon_panel(results)
    b_cells, b_rows, b_cols = build_pyruvate_panel(results)
    c_cells, c_rows, c_cols = build_medium_panel(results)
    c_blocks = _split_blocks(c_cells, c_rows, PANEL_C_BLOCKS)

    margin = 0.18
    gap = 0.60
    a_label_width, b_label_width, c_label_width = 0.80, 1.05, 1.85
    legend_width = 3.10
    header = 0.95

    a_matrix_x = margin + a_label_width
    a_width = len(a_cols) * CELL
    b_matrix_x = a_matrix_x + a_width + gap + b_label_width
    b_width = len(b_cols) * CELL
    legend_x = b_matrix_x + b_width + gap + 0.30

    top_y = header + 0.28
    top_height = max(len(a_rows), len(b_rows)) * CELL
    block_width = len(c_cols) * CELL
    c_matrix_y = top_y + top_height + 0.55 + header
    c_height = max(len(rows) for _, rows in c_blocks) * CELL

    block_content = c_label_width + block_width
    fig_width = (
        max(
            legend_x + legend_width,
            margin + len(c_blocks) * block_content + (len(c_blocks) - 1) * gap,
        )
        + margin
    )
    fig_height = c_matrix_y + c_height + margin + 0.15
    fig = plt.figure(figsize=(fig_width, fig_height))

    ax_a = _axes_inches(fig, a_matrix_x, top_y, a_width, len(a_rows) * CELL)
    _draw_matrix(ax_a, a_cells, a_rows, a_cols, rotate_columns=True)
    _figure_label(fig, margin, top_y - header - 0.24, "A")

    ax_b = _axes_inches(fig, b_matrix_x, top_y, b_width, len(b_rows) * CELL)
    _draw_matrix(ax_b, b_cells, b_rows, b_cols, rotate_columns=True)
    _figure_label(fig, b_matrix_x - b_label_width, top_y - header - 0.24, "B")

    ax_key = _axes_inches(fig, legend_x, top_y - 0.15, 2.30, 0.70)
    _draw_key(ax_key)

    ax_legend = _axes_inches(fig, legend_x, top_y + 0.95, legend_width, 2.05)
    _draw_legend(ax_legend)

    # Spread the blocks across the full width rather than packing them left,
    # so panel C sits under the legend instead of leaving a gutter beside it.
    if len(c_blocks) > 1:
        stride = (fig_width - margin - block_content - margin) / (len(c_blocks) - 1)
    else:
        stride = 0.0

    for index, (block_cells, block_rows) in enumerate(c_blocks):
        block_x = margin + index * stride + c_label_width
        ax_block = _axes_inches(
            fig, block_x, c_matrix_y, block_width, len(block_rows) * CELL
        )
        _draw_matrix(ax_block, block_cells, block_rows, c_cols, rotate_columns=True)
        if index == 0:
            _figure_label(fig, margin, c_matrix_y - header - 0.20, "C")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(
            os.path.join(RESULTS_DIR, f"{output_stem}.{extension}"),
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)

    # The classified table backs the numbers quoted in the caption.
    table = results.drop(columns=["met_id"]).copy()
    table.to_csv(
        os.path.join(RESULTS_DIR, f"{output_stem}_classified.tsv"),
        sep="\t",
        index=False,
    )
    with open(
        os.path.join(RESULTS_DIR, f"{output_stem}_summary.txt"), "w"
    ) as handle:
        handle.write(format_summary(summary) + "\n")
        handle.write(f"Growth threshold: {GROWTH_THRESHOLD} 1/hr\n")

    return results, summary


def main():
    import cobra

    model = cobra.io.read_sbml_model(MODEL_PATH)
    results, _ = generate_phenotype_figure(model)

    discordant = results[results["discordant"]]
    if not discordant.empty:
        print("\nMismatches:")
        print(
            discordant[
                ["minimal_media", "c_source", "growth", "predicted", "category"]
            ].to_string(index=False)
        )

    unrepresentable = results[results["category"] == "no_exchange"]
    if not unrepresentable.empty:
        print("\nNot representable (no exchange reaction):")
        print(
            unrepresentable[["minimal_media", "c_source", "missing_exchanges"]].to_string(
                index=False
            )
        )


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    main()

#!/usr/bin/env python3
"""
Plot Alteromonas MIT1002 secretion rates for nitrogen-containing metabolites
across the diel cycle, one line per f value.

Identifies every exchange reaction in the simulation outputs whose metabolite
formula contains N and whose flux is positive (= secretion under COBRA's
exchange convention) at any (interval, f) combination. Plots one panel per
metabolite.

Use this to test the hypothesis that Amac recycles Pro's nitrogen waste — if
N-containing carbon substrates (glutamate, etc.) are catabolised with biomass
C:N higher than substrate C:N, the model should predict NH3 secretion.
"""

from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results"
FIG_DIR = SCRIPT_DIR / "figs"
FIG_DIR.mkdir(exist_ok=True)

MODEL_FILE = SCRIPT_DIR / "../../model.xml"
FLUX_FILE = RESULTS_DIR / "fluxes_long.csv"

DARK_PERIODS = [(10, 22), (34, 46)]
SECRETION_THRESHOLD = 1e-9  # mmol/gDW/hr; below this, treat as numerical zero

PALETTE = {
    0.0: "#cccccc",
    0.1: "#a6cee3",
    0.5: "#1f78b4",
    1.0: "#33a02c",
    2.0: "#fb9a99",
    5.0: "#e31a1c",
    10.0: "#ff7f00",
}


def shade_dark(ax, alpha: float = 0.15) -> None:
    for d_start, d_end in DARK_PERIODS:
        ax.axvspan(d_start, d_end, color="gray", alpha=alpha, zorder=0, label="_dark")


def contains_n(formula: str | None) -> bool:
    """True if the chemical formula has an N atom (and is not just N-something
    spurious like 'Na' or 'Ni'). Walks character-by-character to be safe."""
    if not formula:
        return False
    i = 0
    while i < len(formula):
        c = formula[i]
        if c == "N":
            # Look ahead: if the next char is lowercase, it's a different element
            # (Na, Ni, Nb, Nd, Ne, etc.)
            if i + 1 < len(formula) and formula[i + 1].islower():
                i += 2
                continue
            return True
        i += 1
    return False


def main() -> None:
    print("Loading model and flux outputs...")
    model = cobra.io.read_sbml_model(str(MODEL_FILE))
    flux_df = pd.read_csv(FLUX_FILE)

    # All exchange reaction IDs in the simulation results
    ex_rxn_ids = sorted(r for r in flux_df["reaction_id"].unique() if r.startswith("EX_"))
    print(f"  {len(ex_rxn_ids)} exchange reactions in flux output")

    # Build {ex_rxn_id: (metabolite_name, formula)} lookup
    info: dict[str, tuple[str, str]] = {}
    for ex_id in ex_rxn_ids:
        if ex_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(ex_id)
        if not rxn.metabolites:
            continue
        met = next(iter(rxn.metabolites))
        info[ex_id] = (met.name or ex_id, met.formula or "")

    # Filter to N-containing
    n_exchanges = {ex: nm for ex, (nm, f) in info.items() if contains_n(f)}
    print(f"  {len(n_exchanges)} of these have N in the formula")

    # Subset fluxes to these reactions and find ones with secretion anywhere
    ex_fluxes = flux_df[flux_df["reaction_id"].isin(n_exchanges)].copy()
    max_flux = ex_fluxes.groupby("reaction_id")["flux"].max()
    secreted = max_flux[max_flux > SECRETION_THRESHOLD].sort_values(ascending=False)
    print(f"  {len(secreted)} N-containing exchanges have positive (secretion) flux:")
    for rxn_id, mx in secreted.items():
        name, formula = info[rxn_id]
        print(f"    {rxn_id:22s}  {name:40s} ({formula})  max={mx:.4g}")

    if len(secreted) == 0:
        print("No N-containing metabolites are secreted — nothing to plot.")
        return

    # Plot: one panel per metabolite, lines = f values
    rxn_ids_to_plot = list(secreted.index)
    n = len(rxn_ids_to_plot)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.2 * nrows), squeeze=False)

    for ax_idx, rxn_id in enumerate(rxn_ids_to_plot):
        ax = axes[ax_idx // ncols, ax_idx % ncols]
        shade_dark(ax)

        sub = ex_fluxes[ex_fluxes["reaction_id"] == rxn_id].copy()
        sub["midpoint_h"] = (sub["interval_start_h"] + sub["interval_end_h"]) / 2

        for f in sorted(sub["f"].unique()):
            fsub = sub[sub["f"] == f].sort_values("midpoint_h")
            ax.plot(
                fsub["midpoint_h"],
                fsub["flux"],
                "o-",
                color=PALETTE.get(f, "#888"),
                linewidth=1.4,
                markersize=3.5,
                label=f"f = {f}",
                zorder=3,
            )

        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_xlim(0, 46)
        ax.set_xticks(range(0, 47, 8))
        ax.tick_params(labelsize=8)

        name, formula = info[rxn_id]
        # Truncate long names
        short = name if len(name) <= 35 else name[:32] + "…"
        ax.set_title(f"{short}\n{rxn_id}  ({formula})", fontsize=9)
        if ax_idx % ncols == 0:
            ax.set_ylabel("Secretion flux\n(mmol gDW⁻¹ hr⁻¹)", fontsize=9)
        if ax_idx // ncols == nrows - 1:
            ax.set_xlabel("Time (h)", fontsize=9)

    # Hide unused axes
    for ax_idx in range(n, nrows * ncols):
        axes[ax_idx // ncols, ax_idx % ncols].set_visible(False)

    # Single legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles = [h for h, l in zip(handles, labels) if not l.startswith("_")]
    labels = [l for l in labels if not l.startswith("_")]
    dark_patch = plt.Rectangle((0, 0), 1, 1, fc="gray", alpha=0.3, label="Dark period")
    fig.legend(
        handles=handles + [dark_patch],
        labels=labels + ["Dark period"],
        loc="lower center",
        ncol=min(8, len(handles) + 1),
        fontsize=9,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(
        "Alteromonas MIT1002 secretion of N-containing metabolites over the diel cycle",
        fontsize=11,
    )
    plt.tight_layout(rect=(0.0, 0.04, 1, 0.97))

    out_path = FIG_DIR / "fig_n_secretions.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()

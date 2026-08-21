"""
Baseline (no forced dissipation) SMF and PMF budgets across carbon sources.

For each substrate we grow the model on an equal-carbon amount (TOTAL_UPTAKE / n_c)
with pFBA, then decompose how the cell BUILDS and SPENDS each ion-motive force
across the membrane:

  SMF -> the Na+_e0 budget    PMF -> the H+_e0 budget

For every TRANSMEMBRANE reaction touching the extracellular ion:
    rate = flux * stoich_coeff(ion_e0)
    rate > 0  -> ion pumped OUT  -> PRODUCES the motive force
    rate < 0  -> ion brought IN  -> CONSUMES (spends) the motive force

Boundary exchange reactions (EX_*) are EXCLUDED on purpose: this figure shows
only transmembrane bioenergetic flux, not the net acid/base exchange the model
uses to balance protons with the medium. As a result the producer and consumer
bars need not sum to exactly equal heights (they differ by that excluded
medium-exchange term, sizeable only for H+).

The story:
  - NaNQR is the sole SMF producer on every substrate (SMF panel, left bars).
  - The Na+/H+ antiporter (bold, same color in both panels) appears as an SMF
    CONSUMER (Na+ in) and a PMF PRODUCER (H+ out) at equal magnitude -- i.e. it
    launders sodium-motive force into proton-motive force. ~1/5 to ~1/3 of the
    PMF driving ATP synthase is sodium-derived.
"""

import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import cobra
import cobra.flux_analysis
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import the shared plot styles from tools/
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
)
from tools.plot_styles import set_plot_style, summer_colors

matplotlib.rcParams.update(
    {
        "font.size": 11,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,  # editable text in Illustrator
        "ps.fonttype": 42,
    }
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SUBSTRATES = [
    {"name": "glucose", "met_id": "cpd00027", "n_c": 6, "na_symport": False},
    {"name": "glycerol", "met_id": "cpd00100", "n_c": 3, "na_symport": False},
    {"name": "glutamate", "met_id": "cpd00023", "n_c": 5, "na_symport": True},
    {"name": "aspartate", "met_id": "cpd00041", "n_c": 4, "na_symport": True},
    {"name": "alanine", "met_id": "cpd00035", "n_c": 3, "na_symport": True},
    {"name": "glycine", "met_id": "cpd00033", "n_c": 2, "na_symport": True},
    {"name": "lysine", "met_id": "cpd00039", "n_c": 6, "na_symport": True},
]

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]

import sys

sys.path.insert(0, str(REPO_ROOT))

from tools.media import MEDIA  # noqa: E402

OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

TOTAL_UPTAKE = 60
BIOMASS_RXN = "bio1_biomass"

# Reaction -> display category. Reactions not listed (and not exchanges) fall
# into "other"; exchanges (EX_*) are dropped entirely (transmembrane only).
CYTOCHROME_RXNS = {"rxn13688_c0", "rxn14412_c0", "rxn14421_c0", "rxn14422_c0"}
AA_SYMPORT_RXNS = {
    "rxn05298_c0",
    "rxn05215_c0",
    "rxn34493_c0",
    "rxn08661_c0",
    "rxn08854_c0",
}

NANQR = "NaNQR (NADH→Na⁺ pump)"
ANTIPORT = "Na⁺/H⁺ antiporter"
ATPSYN = "ATP synthase (F₁)"
CYTO = "Cytochrome ETC"
AASYM = "AA:Na⁺ symporter"


def category(rxn_id):
    if rxn_id.startswith("EX_"):
        return None  # drop boundary exchanges
    if rxn_id == "ec7211_c0":
        return NANQR
    if rxn_id == "rxn05209_c0":
        return ANTIPORT
    if rxn_id == "rxn08173_c0":
        return ATPSYN
    if rxn_id in CYTOCHROME_RXNS:
        return CYTO
    if rxn_id in AA_SYMPORT_RXNS:
        return AASYM
    return "other"


# Palette (antiporter is the bold accent and shared across both panels)
COLORS = {
    ATPSYN: summer_colors["dark_pink"],
    CYTO: summer_colors["pink"],
    ANTIPORT: summer_colors["yellow"],
    NANQR: summer_colors["teal"],
    AASYM: summer_colors["light_blue"],
    "other (produce)": summer_colors["dark_tan"],
    "other (consume)": summer_colors["dark_tan"],
}

# Bottom-to-top stacking order (antiporter at the base, aligned across panels)
STACK_ORDER = [
    ANTIPORT,
    NANQR,
    ATPSYN,
    CYTO,
    AASYM,
    "other (produce)",
    "other (consume)",
]


def run_substrates(model, minimal_media):
    sols, growth = {}, {}
    for sub in SUBSTRATES:
        ex_id = f"EX_{sub['met_id']}_e0"
        with model:
            media = minimal_media.copy()
            media[ex_id] = TOTAL_UPTAKE / sub["n_c"]
            model.medium = media
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sol = cobra.flux_analysis.pfba(model)
        sols[sub["name"]] = sol
        growth[sub["name"]] = sol.fluxes[BIOMASS_RXN]
        print(
            f"  {sub['name']:<10s} growth={sol.fluxes[BIOMASS_RXN]:.4f}  "
            f"NaNQR={sol.fluxes['ec7211_c0']:.2f}"
        )
    return sols, growth


def ion_budget(model, sols, ion_met):
    """DataFrame indexed by substrate, cols = display category, signed rate.
    Exchanges dropped; minor transmembrane terms split into other produce/consume."""
    rows = {}
    for name, sol in sols.items():
        agg = defaultdict(float)
        for r in model.reactions:
            if ion_met not in r.metabolites:
                continue
            cat = category(r.id)
            if cat is None:
                continue
            rate = sol.fluxes[r.id] * r.metabolites[ion_met]
            if abs(rate) < 1e-9:
                continue
            if cat == "other":
                agg["other (produce)" if rate > 0 else "other (consume)"] += rate
            else:
                agg[cat] += rate
        rows[name] = dict(agg)
    df = pd.DataFrame(rows).T.reindex([s["name"] for s in SUBSTRATES]).fillna(0.0)
    return df


def plot_budget(ax, df, title, ylabel):
    subs = df.index.tolist()
    x = np.arange(len(subs))
    width = 0.38

    cols = [c for c in STACK_ORDER if c in df.columns]
    prod_bottom = np.zeros(len(subs))
    cons_bottom = np.zeros(len(subs))
    legend = {}

    for col in cols:
        vals = df[col].values
        prod = np.clip(vals, 0, None)
        cons = np.clip(vals, None, 0)
        color = COLORS.get(col, "#cccccc")
        label = "other" if col.startswith("other") else col
        if prod.any():
            b = ax.bar(
                x - width / 2,
                prod,
                width,
                bottom=prod_bottom,
                color=color,
                edgecolor="white",
                linewidth=0.4,
            )
            prod_bottom += prod
            legend.setdefault(label, b)
        if cons.any():
            b = ax.bar(
                x + width / 2,
                -cons,
                width,
                bottom=cons_bottom,
                color=color,
                edgecolor="white",
                linewidth=0.4,
                hatch="//",
            )
            cons_bottom += -cons
            legend.setdefault(label, b)

    na_set = {s["name"] for s in SUBSTRATES if s["na_symport"]}
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}†" if s in na_set else s for s in subs])
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="center", fontsize=12, pad=8)
    ax.margins(y=0.02)

    # Legend: reaction colors + a produce/consume key
    handles = list(legend.values())
    labels = list(legend.keys())
    produce_key = mpatches.Patch(facecolor="0.75", edgecolor="white", label="produce")
    consume_key = mpatches.Patch(
        facecolor="0.75", edgecolor="white", hatch="//", label="consume"
    )
    handles += [produce_key, consume_key]
    labels += ["produce (left)", "consume (right)"]
    ax.legend(
        handles,
        labels,
        fontsize=8.5,
        frameon=False,
        bbox_to_anchor=(1.01, 1.0),
        loc="upper left",
    )

    # Style
    set_plot_style(ax)


def main():
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
    minimal_media = MEDIA["minimal"]

    print(
        "Running baseline pFBA per substrate (equal carbon, no forced dissipation)..."
    )
    sols, growth = run_substrates(model, minimal_media)

    na_df = ion_budget(model, sols, model.metabolites.cpd00971_e0)
    h_df = ion_budget(model, sols, model.metabolites.cpd00067_e0)

    na_df.to_csv(OUT_PATH / "baseline_smf_budget.csv")
    h_df.to_csv(OUT_PATH / "baseline_pmf_budget.csv")
    pd.Series(growth, name="growth_rate").to_csv(OUT_PATH / "baseline_growth.csv")

    # report the headline: antiporter share of PMF production
    print("\nAntiporter share of PMF production:")
    for s in [x["name"] for x in SUBSTRATES]:
        prod = h_df.loc[s].clip(lower=0).sum()
        share = h_df.loc[s, ANTIPORT] / prod if prod else float("nan")
        print(f"  {s:<10s} {share*100:4.1f}%")

    fig, (ax_smf, ax_pmf) = plt.subplots(2, 1, figsize=(9.5, 9))
    plot_budget(
        ax_smf,
        na_df,
        "Sodium-motive force budget (Na⁺, transmembrane)",
        "Na⁺ flux (mmol gDW⁻¹ h⁻¹)",
    )
    plot_budget(
        ax_pmf,
        h_df,
        "Proton-motive force budget (H⁺, transmembrane)",
        "H⁺ flux (mmol gDW⁻¹ h⁻¹)",
    )
    fig.text(0.01, 0.005, "† Na⁺-symporter substrate", fontsize=8.5, style="italic")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(OUT_PATH / "baseline_smf_pmf_budget.png", dpi=300, bbox_inches="tight")
    print("\nSaved: baseline_smf_pmf_budget.png (300 dpi)")
    print(
        "Saved: baseline_smf_budget.csv, baseline_pmf_budget.csv, baseline_growth.csv"
    )


if __name__ == "__main__":
    main()

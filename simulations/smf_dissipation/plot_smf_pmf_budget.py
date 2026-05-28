"""
Decompose the SMF (Na+) and PMF (H+) budgets across the dissipation sweep
for a single substrate.

For each extracellular ion (Na+_e0, H+_e0), every reaction's contribution to
the e0 pool is rate = flux * stoich_coeff(ion_e0):
  rate > 0  -> reaction pumps the ion OUT of the cell  (PRODUCES the ion-motive force)
  rate < 0  -> reaction brings the ion INTO the cell   (CONSUMES / dissipates the force)

At steady state production == consumption (mass balance), so at each dissipation
level we draw two side-by-side stacked bars: producers (left) and consumers (right),
each segment a reaction. Watching these across the sweep shows the SMF->PMF handoff:
the Na+/H+ antiporter flips from a PMF *producer* (forward, using SMF) to a PMF
*consumer* (reverse, spending PMF to re-export the forced Na+), and ATP synthase's
share of PMF shrinks accordingly.

Modelled on simulations/pro-top10/pmf_vs_smf/plot_smf_pmf_for_single_and_cocktail.py
but (a) sweeps dissipation for one substrate and (b) splits production/consumption
into separate bars instead of +/- halves of one bar.
"""

import warnings
from pathlib import Path

import cobra
import cobra.flux_analysis
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle as pkl

matplotlib.rcParams.update({"font.size": 9, "axes.linewidth": 0.8})

# ---------------------------------------------------------------------------
# Config -- change SUBSTRATE to run for a different carbon source
# ---------------------------------------------------------------------------
SUBSTRATE = {"name": "glutamate", "met_id": "cpd00023", "n_c": 5}
SWEEP_LEVELS = list(range(0, 220, 20))   # dissipation lower bounds to show
RATE_THRESHOLD = 0.4                      # mmol/gDW/hr; smaller contributions -> "other"

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

TOTAL_UPTAKE = 60
FLAGELLA_RXN = "flagella_Na_import_c0"

# Readable labels for the reactions we expect to dominate the ion budgets
LABELS = {
    "ec7211_c0":            "NaNQR (NADH→Na⁺ pump)",
    "rxn05209_c0":          "Na⁺/H⁺ antiporter",
    "rxn05298_c0":          "Glu:Na⁺ symporter",
    "rxn05215_c0":          "Ala:Na⁺ symporter",
    "rxn34493_c0":          "Asp:Na⁺ symporter",
    "rxn08661_c0":          "Gly:Na⁺ symporter",
    "rxn08854_c0":          "Lys:Na⁺ symporter",
    "rxn05313_c0":          "Pi:3Na⁺ symporter",
    FLAGELLA_RXN:           "Flagella Na⁺ import",
    "EX_cpd00971_e0":       "Na⁺ medium exchange",
    "rxn08173_c0":          "ATP synthase (F₁)",
    "rxn13688_c0":          "cyt-c oxidase",
    "rxn14412_c0":          "cyt-c reductase",
    "rxn14421_c0":          "cyt-bo reductase",
    "rxn14422_c0":          "cyt-bo oxidase",
    "EX_cpd00067_e0":       "H⁺ medium exchange",
    "rxn05488_c0":          "acetate:H⁺ symport",
    "rxn05559_c0":          "formate:H⁺ symport",
    "rxn05625_c0":          "nitrite:H⁺ symport",
    "rxn05312_c0":          "Pi:H⁺ export",
    "rxn05209_c0_rev":      "Na⁺/H⁺ antiporter",
}


# Semantic colors for the key players so the antiporter (the SMF<->PMF switcher)
# and the two pumps stand out. Anything not listed falls back to tab20.
COLOR_MAP = {
    "ec7211_c0":      "#e6550d",   # NaNQR  - strong orange (the SMF generator)
    "rxn05209_c0":    "#6a51a3",   # Na+/H+ antiporter - purple (the switcher)
    FLAGELLA_RXN:     "#252525",   # flagella - near black (the imposed load)
    "rxn08173_c0":    "#cb181d",   # ATP synthase - red (PMF -> ATP)
    "rxn05298_c0":    "#31a354",   # Glu symporter - green
    "rxn05215_c0":    "#31a354",   # Ala symporter - green
    "rxn34493_c0":    "#31a354",   # Asp symporter - green
    "rxn08661_c0":    "#31a354",   # Gly symporter - green
    "rxn08854_c0":    "#31a354",   # Lys symporter - green
    "rxn05313_c0":    "#969696",   # Pi:3Na+ symporter - gray
    "rxn13688_c0":    "#3182bd",   # cyt-c oxidase  - blue
    "rxn14412_c0":    "#6baed6",   # cyt-c reductase - light blue
    "rxn14421_c0":    "#08519c",   # cyt-bo reductase - dark blue
    "rxn14422_c0":    "#9ecae1",   # cyt-bo oxidase - pale blue
    "EX_cpd00971_e0": "#d9d9d9",   # Na+ medium exchange - pale gray
    "EX_cpd00067_e0": "#d9d9d9",   # H+ medium exchange - pale gray
}


def short_label(rxn_id: str, model: cobra.Model) -> str:
    if rxn_id in LABELS:
        return LABELS[rxn_id]
    name = model.reactions.get_by_id(rxn_id).name
    return (name[:32] + "…") if len(name) > 33 else name


def build_model() -> cobra.Model:
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
    rxn = cobra.Reaction(FLAGELLA_RXN)
    rxn.name = "Flagella Na+ import (SMF dissipation)"
    rxn.add_metabolites({
        model.metabolites.cpd00971_e0: -1,
        model.metabolites.cpd00971_c0:  1,
    })
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    model.add_reactions([rxn])
    return model


def run_sweep(model, minimal_media):
    """Return {lb: cobra.Solution} for the chosen substrate."""
    ex_id = f"EX_{SUBSTRATE['met_id']}_e0"
    sols = {}
    for lb in SWEEP_LEVELS:
        with model:
            media = minimal_media.copy()
            media[ex_id] = TOTAL_UPTAKE / SUBSTRATE["n_c"]
            model.medium = media
            model.reactions.get_by_id(FLAGELLA_RXN).lower_bound = lb
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sol = cobra.flux_analysis.pfba(model)
            if sol.status == "optimal":
                sols[lb] = sol
            else:
                print(f"  lb={lb}: infeasible, stopping sweep")
                break
    return sols


def ion_budget(model, sols, ion_met):
    """
    Build a DataFrame indexed by dissipation lb, columns = reaction ids,
    values = signed rate of ion appearing in e0 (flux * coeff).
    Positive = produces the force (ion out), negative = consumes it (ion in).
    """
    ion_rxns = [r.id for r in model.reactions if ion_met in r.metabolites]
    rows = {}
    for lb, sol in sols.items():
        row = {}
        for rid in ion_rxns:
            rxn = model.reactions.get_by_id(rid)
            coeff = rxn.metabolites[ion_met]
            rate = sol.fluxes[rid] * coeff
            if abs(rate) > 1e-9:
                row[rid] = rate
        rows[lb] = row
    df = pd.DataFrame(rows).T.fillna(0.0)   # index=lb, cols=rxn
    df = df.loc[:, (df.abs() > 1e-9).any(axis=0)]
    return df


def collapse_small(df, threshold):
    """Move reactions whose max |rate| < threshold into an 'other' column,
    keeping production-other and consumption-other separate."""
    keep_cols, small_cols = [], []
    for c in df.columns:
        if df[c].abs().max() >= threshold:
            keep_cols.append(c)
        else:
            small_cols.append(c)
    out = df[keep_cols].copy()
    if small_cols:
        small = df[small_cols]
        out["other (produce)"] = small.clip(lower=0).sum(axis=1)
        out["other (consume)"] = small.clip(upper=0).sum(axis=1)
    return out


def plot_budget(ax, df, model, title, ylabel):
    """Grouped stacked bars: producers (left) and consumers (right) per lb."""
    lbs = df.index.tolist()
    x = np.arange(len(lbs))
    width = 0.38

    # Order columns by total absolute contribution for stable stacking/colors
    order = df.abs().sum(axis=0).sort_values(ascending=False).index.tolist()
    cmap = plt.get_cmap("tab20")
    fallback_i = 0
    colors = {}
    for col in order:
        if col in COLOR_MAP:
            colors[col] = COLOR_MAP[col]
        elif col.startswith("other"):
            colors[col] = "#bdbdbd"
        else:
            colors[col] = cmap(fallback_i % 20)
            fallback_i += 1

    prod_bottom = np.zeros(len(lbs))
    cons_bottom = np.zeros(len(lbs))
    legend_handles = {}

    for col in order:
        vals = df[col].values
        prod = np.clip(vals, 0, None)     # production (ion out)
        cons = np.clip(vals, None, 0)     # consumption (ion in), negative
        label = "other" if col.startswith("other") else short_label(col, model)

        if prod.any():
            b = ax.bar(x - width / 2, prod, width, bottom=prod_bottom,
                       color=colors[col], edgecolor="white", linewidth=0.3)
            prod_bottom += prod
            legend_handles.setdefault(label, b)
        if cons.any():
            # plot consumption as positive height on the right bar
            b = ax.bar(x + width / 2, -cons, width, bottom=cons_bottom,
                       color=colors[col], edgecolor="white", linewidth=0.3,
                       hatch="///")
            cons_bottom += -cons
            legend_handles.setdefault(label, b)

    ax.set_xticks(x)
    ax.set_xticklabels(lbs)
    ax.set_xlabel("SMF dissipation — forced Na⁺ import (mmol/gDW/hr)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(legend_handles.values(), legend_handles.keys(),
              fontsize=6.5, bbox_to_anchor=(1.01, 1), loc="upper left")
    return df


def main():
    model = build_model()
    with open(REPO_ROOT / "test" / "test_files" / "media" / "media_definitions.pkl", "rb") as fh:
        minimal_media = pkl.load(fh)["minimal"]

    print(f"Running sweep for {SUBSTRATE['name']}...")
    sols = run_sweep(model, minimal_media)

    na_e = model.metabolites.cpd00971_e0
    h_e = model.metabolites.cpd00067_e0

    na_df = collapse_small(ion_budget(model, sols, na_e), RATE_THRESHOLD)
    h_df = collapse_small(ion_budget(model, sols, h_e), RATE_THRESHOLD)

    # Save decompositions (with readable column names)
    na_named = na_df.rename(columns={c: short_label(c, model) for c in na_df.columns
                                     if not c.startswith("other")})
    h_named = h_df.rename(columns={c: short_label(c, model) for c in h_df.columns
                                   if not c.startswith("other")})
    na_named.to_csv(OUT_PATH / f"smf_budget_{SUBSTRATE['name']}.csv")
    h_named.to_csv(OUT_PATH / f"pmf_budget_{SUBSTRATE['name']}.csv")

    # Two-panel figure: SMF (top), PMF (bottom)
    fig, (ax_smf, ax_pmf) = plt.subplots(2, 1, figsize=(11, 10))
    plot_budget(ax_smf, na_df, model,
                f"SMF budget (Na⁺_e0) on {SUBSTRATE['name']}\n"
                "left bar = produces SMF (Na⁺ pumped out) | right bar (hatched) = consumes SMF (Na⁺ in)",
                "Na⁺ flux across membrane (mmol/gDW/hr)")
    plot_budget(ax_pmf, h_df, model,
                f"PMF budget (H⁺_e0) on {SUBSTRATE['name']}\n"
                "left bar = produces PMF (H⁺ pumped out) | right bar (hatched) = consumes PMF (H⁺ in)",
                "H⁺ flux across membrane (mmol/gDW/hr)")
    fig.tight_layout()
    fig.savefig(OUT_PATH / f"smf_pmf_budget_{SUBSTRATE['name']}.pdf", bbox_inches="tight")
    fig.savefig(OUT_PATH / f"smf_pmf_budget_{SUBSTRATE['name']}.png", dpi=150, bbox_inches="tight")
    print(f"Saved: smf_pmf_budget_{SUBSTRATE['name']}.pdf/.png")
    print(f"Saved: smf_budget_{SUBSTRATE['name']}.csv, pmf_budget_{SUBSTRATE['name']}.csv")


if __name__ == "__main__":
    main()

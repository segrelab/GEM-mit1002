"""
Baseline (no forced dissipation) SMF and PMF budgets across carbon sources.

For each substrate we grow the model on an equal-carbon amount (TOTAL_UPTAKE / n_c)
with pFBA, then decompose how the cell BUILDS and SPENDS each ion-motive force:

  SMF  -> the Na+_e0 budget   PMF -> the H+_e0 budget

For every reaction touching the extracellular ion:
    rate = flux * stoich_coeff(ion_e0)
    rate > 0  -> ion pumped OUT  -> PRODUCES the motive force
    rate < 0  -> ion brought IN  -> CONSUMES (spends) the motive force

At steady state production == consumption, so per substrate we draw two
side-by-side stacked bars: producers (left, solid) and consumers (right, hatched).

The point of the figure: Na+-coupled-substrate growth (glutamate, aspartate,
alanine, glycine, lysine) imposes a standing SMF *demand* (the symporters spend
Na+ to import carbon) that NaNQR must continuously regenerate, whereas glucose
and glycerol carry no such Na+ symport demand. This is what motivates having
NaNQR / an explicit Na+ cycle in the model at all -- no unphysiological forcing
required.

Companion to plot_smf_pmf_budget.py (which sweeps forced dissipation for one
substrate); this one compares substrates at their own optimum.
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
# Config
# ---------------------------------------------------------------------------
# Na_symporter flags which substrates ride in on a Na+ symporter (cost the SMF)
SUBSTRATES = [
    {"name": "glucose",   "met_id": "cpd00027", "n_c": 6, "na_symport": False},
    {"name": "glycerol",  "met_id": "cpd00100", "n_c": 3, "na_symport": False},
    {"name": "glutamate", "met_id": "cpd00023", "n_c": 5, "na_symport": True},
    {"name": "aspartate", "met_id": "cpd00041", "n_c": 4, "na_symport": True},
    {"name": "alanine",   "met_id": "cpd00035", "n_c": 3, "na_symport": True},
    {"name": "glycine",   "met_id": "cpd00033", "n_c": 2, "na_symport": True},
    {"name": "lysine",    "met_id": "cpd00039", "n_c": 6, "na_symport": True},
]
RATE_THRESHOLD = 0.4   # mmol/gDW/hr; smaller contributions -> "other"

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

TOTAL_UPTAKE = 60
BIOMASS_RXN = "bio1_biomass"

LABELS = {
    "ec7211_c0":      "NaNQR (NADH→Na⁺ pump)",
    "rxn05209_c0":    "Na⁺/H⁺ antiporter",
    "rxn05298_c0":    "Glu:Na⁺ symporter",
    "rxn05215_c0":    "Ala:Na⁺ symporter",
    "rxn34493_c0":    "Asp:Na⁺ symporter",
    "rxn08661_c0":    "Gly:Na⁺ symporter",
    "rxn08854_c0":    "Lys:Na⁺ symporter",
    "rxn05313_c0":    "Pi:3Na⁺ symporter",
    "EX_cpd00971_e0": "Na⁺ medium exchange",
    "rxn08173_c0":    "ATP synthase (F₁)",
    "rxn13688_c0":    "cyt-c oxidase",
    "rxn14412_c0":    "cyt-c reductase",
    "rxn14421_c0":    "cyt-bo reductase",
    "rxn14422_c0":    "cyt-bo oxidase",
    "EX_cpd00067_e0": "H⁺ medium exchange",
    "rxn05488_c0":    "acetate:H⁺ symport",
    "rxn05559_c0":    "formate:H⁺ symport",
    "rxn05625_c0":    "nitrite:H⁺ symport",
    "rxn05312_c0":    "Pi:H⁺ export",
}

COLOR_MAP = {
    "ec7211_c0":      "#e6550d",   # NaNQR - orange (the SMF generator)
    "rxn05209_c0":    "#6a51a3",   # Na+/H+ antiporter - purple (the switcher)
    "rxn08173_c0":    "#cb181d",   # ATP synthase - red (PMF -> ATP)
    "rxn05298_c0":    "#31a354",   # symporters - green
    "rxn05215_c0":    "#31a354",
    "rxn34493_c0":    "#31a354",
    "rxn08661_c0":    "#31a354",
    "rxn08854_c0":    "#31a354",
    "rxn05313_c0":    "#969696",   # Pi:3Na+ symporter - gray
    "rxn13688_c0":    "#3182bd",   # cyt-c oxidase - blue
    "rxn14412_c0":    "#6baed6",   # cyt-c reductase - light blue
    "rxn14421_c0":    "#08519c",   # cyt-bo reductase - dark blue
    "rxn14422_c0":    "#9ecae1",   # cyt-bo oxidase - pale blue
    "EX_cpd00971_e0": "#d9d9d9",   # Na+ medium exchange - pale gray
    "EX_cpd00067_e0": "#d9d9d9",   # H+ medium exchange - pale gray
}


def short_label(rxn_id, model):
    if rxn_id in LABELS:
        return LABELS[rxn_id]
    name = model.reactions.get_by_id(rxn_id).name
    return (name[:32] + "…") if len(name) > 33 else name


def run_substrates(model, minimal_media):
    """Return {substrate_name: cobra.Solution} at baseline (no forcing)."""
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
        print(f"  {sub['name']:<10s} growth={sol.fluxes[BIOMASS_RXN]:.4f}  "
              f"NaNQR={sol.fluxes['ec7211_c0']:.2f}")
    return sols, growth


def ion_budget(model, sols, ion_met):
    """DataFrame indexed by substrate, cols=reaction ids, signed rate (flux*coeff)."""
    ion_rxns = [r.id for r in model.reactions if ion_met in r.metabolites]
    rows = {}
    for name, sol in sols.items():
        row = {}
        for rid in ion_rxns:
            coeff = model.reactions.get_by_id(rid).metabolites[ion_met]
            rate = sol.fluxes[rid] * coeff
            if abs(rate) > 1e-9:
                row[rid] = rate
        rows[name] = row
    # preserve substrate order
    df = pd.DataFrame(rows).T.reindex([s["name"] for s in SUBSTRATES]).fillna(0.0)
    df = df.loc[:, (df.abs() > 1e-9).any(axis=0)]
    return df


def collapse_small(df, threshold):
    keep, small = [], []
    for c in df.columns:
        (keep if df[c].abs().max() >= threshold else small).append(c)
    out = df[keep].copy()
    if small:
        s = df[small]
        out["other (produce)"] = s.clip(lower=0).sum(axis=1)
        out["other (consume)"] = s.clip(upper=0).sum(axis=1)
    return out


def plot_budget(ax, df, model, title, ylabel):
    subs = df.index.tolist()
    x = np.arange(len(subs))
    width = 0.38

    order = df.abs().sum(axis=0).sort_values(ascending=False).index.tolist()
    cmap = plt.get_cmap("tab20")
    fb = 0
    colors = {}
    for col in order:
        if col in COLOR_MAP:
            colors[col] = COLOR_MAP[col]
        elif col.startswith("other"):
            colors[col] = "#bdbdbd"
        else:
            colors[col] = cmap(fb % 20)
            fb += 1

    prod_bottom = np.zeros(len(subs))
    cons_bottom = np.zeros(len(subs))
    legend = {}

    for col in order:
        vals = df[col].values
        prod = np.clip(vals, 0, None)
        cons = np.clip(vals, None, 0)
        label = "other" if col.startswith("other") else short_label(col, model)
        if prod.any():
            b = ax.bar(x - width / 2, prod, width, bottom=prod_bottom,
                       color=colors[col], edgecolor="white", linewidth=0.3)
            prod_bottom += prod
            legend.setdefault(label, b)
        if cons.any():
            b = ax.bar(x + width / 2, -cons, width, bottom=cons_bottom,
                       color=colors[col], edgecolor="white", linewidth=0.3,
                       hatch="///")
            cons_bottom += -cons
            legend.setdefault(label, b)

    ax.set_xticks(x)
    # mark Na+-symport substrates with a dagger
    na_set = {s["name"] for s in SUBSTRATES if s["na_symport"]}
    ax.set_xticklabels([f"{s}†" if s in na_set else s for s in subs])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(legend.values(), legend.keys(),
              fontsize=6.5, bbox_to_anchor=(1.01, 1), loc="upper left")


def main():
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")
    with open(REPO_ROOT / "test" / "test_files" / "media" / "media_definitions.pkl", "rb") as fh:
        minimal_media = pkl.load(fh)["minimal"]

    print("Running baseline pFBA per substrate (equal carbon, no forced dissipation)...")
    sols, growth = run_substrates(model, minimal_media)

    na_e = model.metabolites.cpd00971_e0
    h_e = model.metabolites.cpd00067_e0

    na_df = collapse_small(ion_budget(model, sols, na_e), RATE_THRESHOLD)
    h_df = collapse_small(ion_budget(model, sols, h_e), RATE_THRESHOLD)

    na_named = na_df.rename(columns={c: short_label(c, model) for c in na_df.columns
                                     if not c.startswith("other")})
    h_named = h_df.rename(columns={c: short_label(c, model) for c in h_df.columns
                                   if not c.startswith("other")})
    na_named.to_csv(OUT_PATH / "baseline_smf_budget.csv")
    h_named.to_csv(OUT_PATH / "baseline_pmf_budget.csv")
    pd.Series(growth, name="growth_rate").to_csv(OUT_PATH / "baseline_growth.csv")

    fig, (ax_smf, ax_pmf) = plt.subplots(2, 1, figsize=(11, 10))
    plot_budget(ax_smf, na_df, model,
                "SMF budget (Na⁺_e0) across carbon sources  (†=Na⁺-symport substrate)\n"
                "left bar = produces SMF (Na⁺ out) | right bar (hatched) = consumes SMF (Na⁺ in)",
                "Na⁺ flux across membrane (mmol/gDW/hr)")
    plot_budget(ax_pmf, h_df, model,
                "PMF budget (H⁺_e0) across carbon sources\n"
                "left bar = produces PMF (H⁺ out) | right bar (hatched) = consumes PMF (H⁺ in)",
                "H⁺ flux across membrane (mmol/gDW/hr)")
    fig.tight_layout()
    fig.savefig(OUT_PATH / "baseline_smf_pmf_budget.pdf", bbox_inches="tight")
    fig.savefig(OUT_PATH / "baseline_smf_pmf_budget.png", dpi=150, bbox_inches="tight")
    print("Saved: baseline_smf_pmf_budget.pdf/.png")
    print("Saved: baseline_smf_budget.csv, baseline_pmf_budget.csv, baseline_growth.csv")


if __name__ == "__main__":
    main()

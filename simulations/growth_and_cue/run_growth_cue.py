"""MIT1002 growth rate and carbon-use efficiency (CUE) across a broad substrate panel.

Runs pFBA on a curated single-substrate panel at fixed total carbon uptake
(60 mmol C / gDW / hr) and reports, per substrate:
  - growth rate (biomass flux)
  - CUE = 1 - CO2_secreted / C_taken_up

Bars are coloured by the substrate's entry point into central metabolism and
sorted by growth rate. Self-contained: re-runs the simulations, no dependence
on the (retired) PCA analysis.
"""

import re
import warnings
from pathlib import Path
from typing import Optional
import pickle as pkl
import sys

import cobra
import cobra.flux_analysis
from gem_utilities import media as media_utils
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
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
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

# Import plot_styles.py from the root of the repo
sys.path.append(str(REPO_ROOT))
from plot_styles import summer_colors

TOTAL_UPTAKE = 60  # mmol C / gDW / hr
BIOMASS_RXN = "bio1_biomass"
CO2_EX_RXN = "EX_cpd00011_e0"

# Exchange metabolites whose max |flux| across substrates is below this are
# lumped into a single grey "Other" segment (keeps the trace-ion colours out).
EX_FLUX_THRESHOLD = 1.0  # mmol / gDW / hr

# Color palette for exchange metabolites
# The "Summer" color palette with a few extra colors to avoid repeats
# Remove colors that I hard set below
EX_PALETTE = [
    v
    for k, v in summer_colors.items()
    if k not in {"yellow", "light_blue", "dark_tan", "light_tan"}
]
EX_PALETTE += [
    "#9B6A9E",  # mauve
    "#3F6F8C",  # steel blue
    "#C98B8B",  # dusty rose
    "#8C4F6B",  # plum
    "#6E7B8B",  # slate
    "#BFB13A",  # brass
]
# Define an eye-catching color for ammonium, since that is the primary result
NH3_COLOR = summer_colors["yellow"]
# Define muted colors for metabolites that are unimportant
H_COLOR = summer_colors["light_blue"]
CO2_COLOR = summer_colors["dark_tan"]
H2O_COLOR = summer_colors["light_tan"]
# Define a single color for the carbon source
# So there aren't unique colros for each bar
# OK to be the same color as CO2, since the carbon source is ony ever taken up
# And CO2 is only ever released
CARBON_SOURCE_COLOR = summer_colors["dark_tan"]  # all substrate carbon sources
OTHER_COLOR = "#bdbdbd"  # grey — collapsed trace metabolites

# Entry point into central metabolism = the first central-metabolism
# intermediate the substrate's catabolism produces.
ENTRY_POINT = {
    "cpd00027": "Glycolysis",  # Glucose      → G6P
    "cpd00108": "Glycolysis",  # Galactose    → G1P → G6P
    "cpd00080": "Glycolysis",  # Glycerol-3-P → DHAP
    "cpd00035": "Pyruvate",  # Alanine      → Pyr
    "cpd00033": "Pyruvate",  # Glycine      → Ser → Pyr
    "cpd23538": "Pyruvate",  # DHPS         → Pyr + sulfite
    "cpd00029": "Acetyl-CoA",  # Acetate      → AcCoA
    "cpd00797": "Acetyl-CoA",  # 3-HB         → 2× AcCoA
    "cpd00107": "Acetyl-CoA",  # Leucine      → AcCoA
    "cpd00039": "Acetyl-CoA",  # Lysine       → AcCoA
    "cpd00023": "TCA — α-KG",  # Glutamate    → α-KG
    "cpd00129": "TCA — α-KG",  # Proline      → Glu → α-KG
    "cpd00051": "TCA — α-KG",  # Arginine     → Glu → α-KG
    "cpd00041": "TCA — OAA",  # Aspartate    → OAA
    "cpd00156": "TCA — Succinyl-CoA",  # Valine       → succinyl-CoA
    "cpd00322": "TCA — Succinyl-CoA",  # Isoleucine   → succinyl-CoA
    "cpd00123": "TCA — Succinyl-CoA",  # KIC          → succinyl-CoA
    "cpd00069": "Aromatic catabolism",  # Tyrosine     → fumarate + AcAcCoA
    "cpd00127": "Aromatic catabolism",  # Phenol       → succinyl-CoA + AcCoA
}

# SMF-figure palette, extended with two harmonising tones for the 7 categories.
ENTRY_POINT_COLORS = {
    "Glycolysis": summer_colors["light_blue"],
    "Pyruvate": summer_colors["teal"],
    "Acetyl-CoA": summer_colors["yellow"],
    "TCA — α-KG": summer_colors["pink"],
    "TCA — OAA": summer_colors["dark_pink"],
    "TCA — Succinyl-CoA": summer_colors["green"],
    "Aromatic catabolism": "#9B6A9E",  # muted mauve (added)
}
ENTRY_POINT_ORDER = list(ENTRY_POINT_COLORS)


def count_carbons(formula: str) -> Optional[int]:
    """Number of carbon atoms in a molecular formula string."""
    if not formula:
        return None
    m = re.search(r"C(\d*)", formula)
    if m:
        n = m.group(1)
        return int(n) if n else 1
    return 0


def load_substrates(model: cobra.Model, media_defs: dict) -> pd.DataFrame:
    """Build the single-substrate panel from the known-growth-phenotypes TSV."""
    df = pd.read_csv(TEST_FILE_DIR / "known_growth_phenotypes.tsv", sep="\t")
    df = df[df["growth"] == "Yes"].copy()
    df = df[~df["met_id"].astype(str).str.contains(",", na=True)].copy()

    if "pro_exomet" in df.columns:
        df["pro_exomet"] = df["pro_exomet"].fillna("No")
    else:
        df["pro_exomet"] = "No"

    # Keep the Pro-exometabolite-annotated row when a substrate appears in
    # multiple media contexts (stable sort puts "Yes" first, then dedup).
    df = (
        df.sort_values(
            "pro_exomet",
            key=lambda s: s.map({"Yes": 0}).fillna(1),
            ascending=True,
            kind="stable",
        )
        .drop_duplicates(subset="met_id", keep="first")
        .copy()
    )

    # Aspartate and glycine were measured in the Pro exometabolome; add if absent.
    for name, met_id in [("Aspartate", "cpd00041"), ("Glycine", "cpd00033")]:
        if met_id not in df["met_id"].values:
            df = pd.concat(
                [
                    df,
                    pd.DataFrame(
                        [
                            {
                                "minimal_media": "l1",
                                "c_source": name,
                                "met_id": met_id,
                                "growth": "Yes",
                                "pro_exomet": "Yes",
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

    rxn_ids = {r.id for r in model.reactions}
    records = []
    for _, row in df.iterrows():
        met_id = str(row["met_id"]).strip()
        c_source = str(row["c_source"]).strip()
        ex_id = f"EX_{met_id}_e0"
        media_key = str(row["minimal_media"]).strip()

        if ex_id not in rxn_ids:
            print(f"  SKIP (no exchange rxn)  : {c_source} ({met_id})")
            continue
        if media_key not in media_defs:
            print(f"  SKIP (unknown media '{media_key}'): {c_source}")
            continue

        try:
            n_c = count_carbons(model.metabolites.get_by_id(f"{met_id}_e0").formula)
        except KeyError:
            n_c = None
            for met in model.metabolites:
                if met.id.startswith(met_id):
                    n_c = count_carbons(met.formula)
                    break
        if not n_c:
            print(f"  SKIP (can't determine n_c): {c_source}")
            continue

        records.append(
            {
                "name": c_source,
                "met_id": met_id,
                "exchange_id": ex_id,
                "media_key": media_key,
                "n_c": n_c,
                "entry_point": ENTRY_POINT.get(met_id, "Other"),
            }
        )

    substrate_df = pd.DataFrame(records)
    print(f"\nSubstrate panel: {len(substrate_df)} substrates")
    return substrate_df


def run_pfba(model, media_defs, substrate_df) -> tuple:
    """Run pFBA per substrate.

    Returns (summary_df, ex_records) where summary_df holds growth + CUE and
    ex_records maps {substrate: {exchange_rxn_id: flux}} (non-zero fluxes only).
    """
    ex_rxn_ids = [r.id for r in model.reactions if r.id.startswith("EX_")]
    rows = []
    ex_records = {}
    for _, row in substrate_df.iterrows():
        name = row["name"]
        media = media_utils.clean_media(model, media_defs[row["media_key"]])
        media[row["exchange_id"]] = TOTAL_UPTAKE / row["n_c"]
        with model:
            model.medium = media
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sol = cobra.flux_analysis.pfba(model)
                growth = sol.fluxes[BIOMASS_RXN]
                if growth < 1e-6:
                    print(f"  WARN (near-zero growth)  : {name}  mu={growth:.5f}")
                    continue
                co2 = sol.fluxes.get(CO2_EX_RXN, 0.0)
                cue = 1.0 - (co2 / TOTAL_UPTAKE)
                rows.append(
                    {
                        "substrate": name,
                        "met_id": row["met_id"],
                        "growth_rate": growth,
                        "co2_flux": co2,
                        "cue": cue,
                        "entry_point": row["entry_point"],
                    }
                )
                ex_records[name] = {
                    rid: sol.fluxes[rid]
                    for rid in ex_rxn_ids
                    if abs(sol.fluxes[rid]) > 1e-9
                }
                print(f"  OK   {name:28s}  mu={growth:.4f}  CUE={cue:.3f}")
            except Exception as exc:
                print(f"  FAIL {name}: {exc}")
    return pd.DataFrame(rows).set_index("substrate"), ex_records


def plot_growth_cue(summary_df: pd.DataFrame, out_path: Path) -> None:
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
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_growth_cue_correlation(summary_df: pd.DataFrame, out_path: Path) -> None:
    df = summary_df.reset_index()
    colors = [ENTRY_POINT_COLORS.get(c, "#9E9E9E") for c in df["entry_point"]]
    x = df["growth_rate"].values.astype(float)
    y = df["cue"].values.astype(float)

    r, p = stats.pearsonr(x, y)
    rho, p_s = stats.spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.scatter(x, y, c=colors, s=120, edgecolor="white", linewidth=0.6, zorder=3)

    # Least-squares fit line
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, slope * xs + intercept, color="#555555", ls="--", lw=1.2, zorder=2)

    # Substrate labels. Default is offset to the right of the point; a few
    # near-coincident points get a manual override so labels don't overlap.
    label_offsets = {  # substrate: (dx, dy, ha, va)
        "Galactose": (0, 9, "center", "bottom"),  # sits almost on Glucose
    }
    for _, row in df.iterrows():
        dx, dy, ha, va = label_offsets.get(row["substrate"], (5, 0, "left", "center"))
        ax.annotate(
            row["substrate"],
            (row["growth_rate"], row["cue"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7,
            ha=ha,
            va=va,
            color="#555555",
        )

    ax.set_xlabel("Growth rate (h⁻¹)")
    ax.set_ylabel("Carbon-use efficiency")
    ax.set_title(
        "Growth rate vs. carbon-use efficiency across substrates", fontsize=12, pad=8
    )

    # Correlation metrics
    txt = (
        f"Pearson  r = {r:.2f}  (p = {p:.2g})\n"
        f"Spearman ρ = {rho:.2f}  (p = {p_s:.2g})"
    )
    ax.text(
        0.025,
        0.025,
        txt,
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="left",
        linespacing=1.4,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc", alpha=0.9),
    )

    # Legend: entry-point colors + the fit line
    handles = [
        mpatches.Patch(facecolor=ENTRY_POINT_COLORS[lbl], edgecolor="white", label=lbl)
        for lbl in ENTRY_POINT_ORDER
        if lbl in df["entry_point"].values
    ]
    handles.append(
        Line2D([0], [0], color="#555555", ls="--", lw=1.2, label="linear fit")
    )
    ax.legend(
        handles=handles,
        title="Entry point into\ncentral metabolism",
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
        title_fontsize=9.5,
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}  (Pearson r={r:.2f}, p={p:.2g})")


def build_exchange_df(model, ex_records, substrate_order, threshold):
    """{substrate: {ex_rxn: flux}} -> DataFrame (substrate x metabolite name).

    Columns are renamed from exchange-reaction id to the metabolite name.
    Trace metabolites (max |flux| < threshold) are collapsed into 'Other'.
    """
    df = pd.DataFrame(ex_records).T.reindex(substrate_order).fillna(0.0)
    rename = {}
    for rid in df.columns:
        met = next(iter(model.reactions.get_by_id(rid).metabolites))
        rename[rid] = met.name
    df = df.rename(columns=rename)
    # collapse duplicate metabolite-name columns if any arise (version-safe)
    if df.columns.duplicated().any():
        df = df.T.groupby(level=0).sum().T

    if threshold > 0:
        keep = [c for c in df.columns if df[c].abs().max() >= threshold]
        small = [c for c in df.columns if c not in keep]
        out = df[keep].copy()
        if small:
            out["Other"] = df[small].sum(axis=1)
        df = out
    return df


def _merge_carbon_sources(df, carbon_source_names):
    """Collapse all substrate-carbon-source columns into one 'Carbon source'."""
    cs = [c for c in df.columns if c in carbon_source_names]
    rest = [c for c in df.columns if c not in carbon_source_names]
    out = df[rest].copy()
    if cs:
        out["Carbon source"] = df[cs].sum(axis=1)
    return out


def _ordered_cols(df, pin_first=None, pin_last=None):
    """Column order: pinned-first, then regular by descending magnitude, pinned-last."""
    skip = {pin_first, pin_last}
    middle = sorted(
        (c for c in df.columns if c not in skip),
        key=lambda c: df[c].abs().sum(),
        reverse=True,
    )
    head = [pin_first] if pin_first in df.columns else []
    tail = [pin_last] if pin_last in df.columns else []
    return head + middle + tail


def _stacked_bar(df, colors, title, ylabel, out_path):
    fig, ax = plt.subplots(figsize=(13, 7))
    x = np.arange(len(df.index))
    bottom = np.zeros(len(df))
    for col in df.columns:
        vals = df[col].values
        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=colors[col],
            label=col,
            edgecolor="white",
            linewidth=0.3,
        )
        bottom += vals
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12, pad=8)
    ax.margins(x=0.01)
    ax.legend(
        title="Metabolite",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=7,
        title_fontsize=8.5,
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_exchange_stacks(ex_df, carbon_source_names, out_dir):
    """Two stacked-bar charts, both drawn upward, with shared metabolite colours.

    Uptake fluxes are negated so they read as positive bars. All substrate
    carbon sources are merged into one 'Carbon source' segment in the uptake
    chart (each appears on a single bar). Byproducts/co-substrates keep
    consistent colours across both charts.
    """
    uptake = -ex_df.clip(upper=0)  # negative fluxes -> positive uptake
    exud = ex_df.clip(lower=0)  # positive fluxes -> exudation
    uptake = uptake.loc[:, (uptake != 0).any(axis=0)]
    exud = exud.loc[:, (exud != 0).any(axis=0)]
    uptake = _merge_carbon_sources(uptake, carbon_source_names)

    # Drop H2O from the exudation chart
    exud = exud.drop(columns=["H2O [e0]"], errors="ignore")

    # Shared colour map over the "regular" metabolites (everything except the
    # fixed-colour Carbon source / Other), ordered by total magnitude so shared
    # metabolites get the same colour in both charts.
    special = {"H2O [e0]", "H+ [e0]", "CO2 [e0]", "NH3 [e0]", "Carbon source", "Other"}
    regular = (set(uptake.columns) | set(exud.columns)) - special

    def total_mag(c):
        m = 0.0
        if c in uptake.columns:
            m += uptake[c].abs().sum()
        if c in exud.columns:
            m += exud[c].abs().sum()
        return m

    ordered_reg = sorted(regular, key=total_mag, reverse=True)
    colors = {c: EX_PALETTE[i % len(EX_PALETTE)] for i, c in enumerate(ordered_reg)}
    # Define the special colors
    colors["H+ [e0]"] = H_COLOR
    colors["H2O [e0]"] = H2O_COLOR
    colors["CO2 [e0]"] = CO2_COLOR
    colors["NH3 [e0]"] = NH3_COLOR
    colors["Carbon source"] = CARBON_SOURCE_COLOR
    colors["Other"] = OTHER_COLOR

    up_cols = _ordered_cols(uptake, pin_first="Carbon source", pin_last="Other")
    ex_cols = _ordered_cols(exud, pin_last="Other")
    _stacked_bar(
        uptake[up_cols],
        colors,
        "MIT1002 uptake fluxes across substrates",
        "Uptake flux (mmol gDW⁻¹ h⁻¹)",
        out_dir / "uptake_fluxes.png",
    )
    _stacked_bar(
        exud[ex_cols],
        colors,
        "MIT1002 exudation fluxes across substrates",
        "Exudation flux (mmol gDW⁻¹ h⁻¹)",
        out_dir / "exudation_fluxes.png",
    )


def main():
    print("Loading model...")
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    print("Loading media definitions...")
    with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
        media_defs = pkl.load(f)

    print("\nBuilding substrate panel...")
    substrate_df = load_substrates(model, media_defs)

    print("\nRunning pFBA simulations...")
    summary_df, ex_records = run_pfba(model, media_defs, substrate_df)
    print(f"\nSuccessful: {len(summary_df)}/{len(substrate_df)} substrates")

    summary_df.to_csv(OUT_PATH / "growth_and_cue.csv")
    print("\nPlotting growth rate and CUE...")
    plot_growth_cue(summary_df, OUT_PATH / "growth_and_cue.png")
    plot_growth_cue_correlation(summary_df, OUT_PATH / "growth_vs_cue_scatter.png")

    # Exchange (uptake/exudation) stacked bars, substrates ordered by growth rate
    print("\nPlotting exchange fluxes...")
    order = summary_df.sort_values("growth_rate", ascending=False).index.tolist()
    ex_df = build_exchange_df(model, ex_records, order, EX_FLUX_THRESHOLD)
    ex_df.to_csv(OUT_PATH / "exchange_fluxes.csv")
    # Metabolite names that are substrate carbon sources (merged in uptake chart)
    carbon_source_names = {
        next(iter(model.reactions.get_by_id(r["exchange_id"]).metabolites)).name
        for _, r in substrate_df.iterrows()
    }
    plot_exchange_stacks(ex_df, carbon_source_names, OUT_PATH)
    print(f"\nAll results saved to: {OUT_PATH}")


if __name__ == "__main__":
    main()

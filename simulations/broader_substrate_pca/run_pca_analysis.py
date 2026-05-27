"""Broader substrate panel PCA of MIT1002 flux distributions.

Runs pFBA on a curated substrate panel at fixed total carbon uptake
(60 mmol C/gDW/hr), builds a (substrate x reaction) flux matrix, and runs:
  1. Standardized PCA — patterns of relative variation, magnitude-independent.
  2. Growth-rate-normalized then standardized PCA — removes the growth-rate
     axis so PC1 reflects metabolic routing rather than overall activity.
"""

from pathlib import Path
import pickle as pkl
import re
import warnings
from typing import Optional

import cobra
from gem_utilities import media as media_utils
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

TOTAL_UPTAKE = 60       # mmol C / gDW / hr
BIOMASS_RXN = "bio1_biomass"
CO2_EX_RXN = "EX_cpd00011_e0"

# Metabolic entry point into central metabolism, keyed by ModelSEED met_id.
# Entry point = the central-metabolism intermediate the substrate first produces
# (or the pathway step at which it joins glycolysis / TCA cycle).
ENTRY_POINT = {
    # ── Glycolysis ──────────────────────────────────────────────────────────
    "cpd00027": "Glycolysis",        # Glucose       → G6P
    "cpd00108": "Glycolysis",        # Galactose     → G1P → G6P (Leloir pathway)
    "cpd00080": "Glycolysis",        # Glycerol-3-P  → DHAP (aldolase/G3P dehydrogenase)
    # ── Pyruvate ────────────────────────────────────────────────────────────
    "cpd00035": "Pyruvate",          # Alanine       → Pyr (alanine aminotransferase)
    "cpd00033": "Pyruvate",          # Glycine       → Ser → Pyr (serine hydroxymethyltransferase)
    "cpd23538": "Pyruvate",          # DHPS          → C3 sulfo-intermediate → Pyr + sulfite
    # ── Acetyl-CoA ──────────────────────────────────────────────────────────
    "cpd00029": "Acetyl-CoA",        # Acetate       → AcCoA (acetyl-CoA synthetase)
    "cpd00797": "Acetyl-CoA",        # 3-HB          → AcAcCoA → 2× AcCoA
    "cpd00107": "Acetyl-CoA",        # Leucine       → HMG-CoA → AcCoA (ketogenic only)
    "cpd00039": "Acetyl-CoA",        # Lysine        → saccharopine/pipecolate → AcCoA
    # ── TCA — α-ketoglutarate ───────────────────────────────────────────────
    "cpd00023": "TCA — α-KG",        # Glutamate     → α-KG (glutamate dehydrogenase)
    "cpd00129": "TCA — α-KG",        # Proline       → Glu → α-KG
    "cpd00051": "TCA — α-KG",        # Arginine      → Glu → α-KG (arginine succinyltransferase)
    # ── TCA — oxaloacetate ──────────────────────────────────────────────────
    "cpd00041": "TCA — OAA",         # Aspartate     → OAA (aspartate aminotransferase)
    # ── TCA — succinyl-CoA (via propionyl-CoA) ──────────────────────────────
    "cpd00156": "TCA — Succinyl-CoA",  # Valine      → propionyl-CoA → succinyl-CoA
    "cpd00322": "TCA — Succinyl-CoA",  # Isoleucine  → propionyl-CoA → succinyl-CoA (+AcCoA)
    "cpd00123": "TCA — Succinyl-CoA",  # KIC         → isobutyryl-CoA → propionyl-CoA → succinyl-CoA
    # ── Aromatic catabolism (split TCA entry) ───────────────────────────────
    "cpd00069": "Aromatic catabolism", # Tyrosine    → homogentisate → fumarate + AcAcCoA
    "cpd00127": "Aromatic catabolism", # Phenol      → catechol → β-ketoadipate → succinyl-CoA + AcCoA
}

ENTRY_POINT_COLORS = {
    "Glycolysis":           "#1f77b4",  # blue
    "Pyruvate":             "#ff7f0e",  # orange
    "Acetyl-CoA":           "#2ca02c",  # green
    "TCA — α-KG":           "#d62728",  # red
    "TCA — OAA":            "#9467bd",  # purple
    "TCA — Succinyl-CoA":   "#8c564b",  # brown
    "Aromatic catabolism":  "#e377c2",  # pink
}


# ── Helpers ──────────────────────────────────────────────────────────────────

def count_carbons(formula: str) -> Optional[int]:
    """Return number of carbon atoms from a molecular formula string."""
    if not formula:
        return None
    m = re.search(r"C(\d*)", formula)
    if m:
        n = m.group(1)
        return int(n) if n else 1
    return 0


# ── Substrate loading ─────────────────────────────────────────────────────────

def load_substrates(model: cobra.Model, media_defs: dict) -> pd.DataFrame:
    """Build the substrate panel from the known-growth-phenotypes TSV.

    Filters for confirmed single-substrate growth, deduplicates by met_id,
    appends aspartate and glycine, then verifies exchange reactions exist.
    """
    tsv = TEST_FILE_DIR / "known_growth_phenotypes.tsv"
    df = pd.read_csv(tsv, sep="\t")

    # Keep only confirmed-growth, single-carbon-source rows
    df = df[df["growth"] == "Yes"].copy()
    df = df[~df["met_id"].astype(str).str.contains(",", na=True)].copy()

    # First occurrence wins when the same met_id appears in multiple media
    df = df.drop_duplicates(subset="met_id", keep="first").copy()

    # Manually add aspartate and glycine (absent from TSV; use l1 as base)
    for name, met_id in [("Aspartate", "cpd00041"), ("Glycine", "cpd00033")]:
        if met_id not in df["met_id"].values:
            df = pd.concat(
                [df, pd.DataFrame([{"minimal_media": "l1", "c_source": name, "met_id": met_id, "growth": "Yes"}])],
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

        # Carbon count from extracellular metabolite formula
        try:
            met = model.metabolites.get_by_id(f"{met_id}_e0")
            n_c = count_carbons(met.formula)
        except KeyError:
            n_c = None
            for met in model.metabolites:
                if met.id.startswith(met_id):
                    n_c = count_carbons(met.formula)
                    break

        if not n_c:
            print(f"  SKIP (can't determine n_c): {c_source}")
            continue

        records.append({
            "name":        c_source,
            "met_id":      met_id,
            "exchange_id": ex_id,
            "media_key":   media_key,
            "n_c":         n_c,
            "entry_point": ENTRY_POINT.get(met_id, "Other"),
        })

    substrate_df = pd.DataFrame(records)
    print(f"\nSubstrate panel: {len(substrate_df)} substrates")
    print(
        substrate_df[["name", "met_id", "n_c", "media_key", "entry_point"]]
        .to_string(index=False)
    )
    return substrate_df


# ── pFBA simulations ──────────────────────────────────────────────────────────

def run_pfba(
    model: cobra.Model, media_defs: dict, substrate_df: pd.DataFrame
) -> tuple:
    """Run pFBA for each substrate.

    Returns:
        flux_matrix : DataFrame (substrate x reaction), index = substrate name
        summary_df  : DataFrame with growth_rate, co2_flux, cue per substrate
    """
    flux_records = {}
    summary_records = []

    for _, row in substrate_df.iterrows():
        name = row["name"]
        base_media = media_defs[row["media_key"]]
        media = media_utils.clean_media(model, base_media)
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

                co2_flux = sol.fluxes.get(CO2_EX_RXN, 0.0)
                cue = 1.0 - (co2_flux / TOTAL_UPTAKE)

                flux_records[name] = sol.fluxes.to_dict()
                summary_records.append({
                    "substrate":   name,
                    "met_id":      row["met_id"],
                    "growth_rate": growth,
                    "co2_flux":    co2_flux,
                    "cue":         cue,
                    "entry_point": row["entry_point"],
                })
                print(f"  OK   {name:28s}  mu={growth:.4f}  CUE={cue:.3f}")

            except Exception as exc:
                print(f"  FAIL {name}: {exc}")

    flux_matrix = pd.DataFrame(flux_records).T.fillna(0)
    summary_df = pd.DataFrame(summary_records).set_index("substrate")
    return flux_matrix, summary_df


# ── PCA utilities ─────────────────────────────────────────────────────────────

def prepare_matrix(
    flux_matrix: pd.DataFrame,
    growth_rates: Optional[pd.Series] = None,
) -> tuple:
    """Drop zero/constant columns, optionally normalize by growth rate, standardize.

    Returns (scaled_array, reaction_names, substrate_names).
    """
    mat = flux_matrix.copy()

    if growth_rates is not None:
        mat = mat.div(growth_rates, axis=0)

    # Drop reactions that are zero in every condition
    mat = mat.loc[:, (mat != 0).any(axis=0)]

    scaler = StandardScaler()
    scaled = scaler.fit_transform(mat.values)

    # Drop columns that are NaN after scaling (constant across conditions)
    valid = ~np.isnan(scaled).any(axis=0)
    scaled = scaled[:, valid]
    reactions = mat.columns[valid].tolist()

    return scaled, reactions, mat.index.tolist()


def run_pca(
    scaled: np.ndarray, substrates: list, reactions: list, n_components: int = 5
) -> tuple:
    """Fit PCA; return (scores_df, loadings_df, explained_variance_ratio)."""
    n = min(n_components, scaled.shape[0] - 1, scaled.shape[1])
    pca = PCA(n_components=n)
    scores = pca.fit_transform(scaled)
    pcs = [f"PC{i+1}" for i in range(n)]
    scores_df = pd.DataFrame(scores, index=substrates, columns=pcs)
    loadings_df = pd.DataFrame(pca.components_.T, index=reactions, columns=pcs)
    return scores_df, loadings_df, pca.explained_variance_ratio_


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_scores(
    scores_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    explained_var: np.ndarray,
    title: str,
    out_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 8))

    classes = summary_df.loc[scores_df.index, "entry_point"]
    for cls in sorted(classes.unique()):
        subs = scores_df[classes == cls]
        ax.scatter(
            subs["PC1"], subs["PC2"],
            color=ENTRY_POINT_COLORS.get(cls, "#9E9E9E"),
            s=180, edgecolor="black", linewidth=0.6,
            label=cls, zorder=3,
        )
        for sub, row in subs.iterrows():
            ax.annotate(
                sub, (row["PC1"], row["PC2"]),
                xytext=(7, 0), textcoords="offset points",
                fontsize=8, va="center", color="#444444",
            )

    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.axvline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel(f"PC1 ({explained_var[0]:.1%} variance)", fontsize=12)
    ax.set_ylabel(f"PC2 ({explained_var[1]:.1%} variance)", fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(title="Class", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_loadings(
    loadings_df: pd.DataFrame,
    model: cobra.Model,
    explained_var: np.ndarray,
    title_prefix: str,
    out_path: Path,
    top_n: int = 15,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16, max(6, 0.4 * top_n)))

    for i, pc in enumerate(["PC1", "PC2"]):
        top_idx = loadings_df[pc].abs().nlargest(top_n).index
        loadings = loadings_df.loc[top_idx, pc].sort_values()

        labels = []
        for rxn_id in loadings.index:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                raw = rxn.name or rxn_id
                label = (raw[:38] + "…") if len(raw) > 38 else raw
            except KeyError:
                label = rxn_id
            labels.append(label)

        colors = ["#d95f0e" if x < 0 else "#2c7fb8" for x in loadings.values]
        axes[i].barh(range(len(loadings)), loadings.values, color=colors)
        axes[i].set_yticks(range(len(loadings)))
        axes[i].set_yticklabels(labels, fontsize=8)
        axes[i].axvline(0, color="black", lw=0.5)
        axes[i].set_xlabel(f"{pc} loading", fontsize=11)
        axes[i].set_title(
            f"{title_prefix} — top {top_n} by |loading| on {pc}\n"
            f"({explained_var[i]:.1%} variance explained)",
            fontsize=10,
        )
        sns.despine(ax=axes[i])

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_growth_cue(summary_df: pd.DataFrame, out_path: Path) -> None:
    df = summary_df.sort_values("growth_rate", ascending=False).reset_index()
    colors = [ENTRY_POINT_COLORS.get(c, "#9E9E9E") for c in df["entry_point"]]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    axes[0].bar(df["substrate"], df["growth_rate"], color=colors, edgecolor="black", lw=0.5)
    axes[0].set_ylabel("Growth rate (1/hr)", fontsize=11)
    axes[0].set_title("MIT1002 growth rate and CUE — broader substrate panel", fontsize=12)
    sns.despine(ax=axes[0])

    axes[1].bar(df["substrate"], df["cue"], color=colors, edgecolor="black", lw=0.5)
    axes[1].set_ylabel("CUE", fontsize=11)
    axes[1].set_ylim(0, 1)
    plt.setp(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=9)
    sns.despine(ax=axes[1])

    handles = [
        mpatches.Patch(color=c, label=lbl)
        for lbl, c in ENTRY_POINT_COLORS.items()
        if lbl in df["entry_point"].values
    ]
    axes[0].legend(
        handles=handles, title="Class",
        bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9,
    )

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading model...")
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    print("Loading media definitions...")
    with open(TEST_FILE_DIR / "media" / "media_definitions.pkl", "rb") as f:
        media_defs = pkl.load(f)

    print("\nBuilding substrate panel...")
    substrate_df = load_substrates(model, media_defs)

    print("\nRunning pFBA simulations...")
    flux_matrix, summary_df = run_pfba(model, media_defs, substrate_df)

    n_ok = len(summary_df)
    n_total = len(substrate_df)
    print(f"\nSuccessful: {n_ok}/{n_total} substrates")

    if n_ok < 2:
        raise RuntimeError("Too few successful simulations to run PCA.")

    flux_matrix.to_csv(OUT_PATH / "flux_matrix.csv")
    summary_df.to_csv(OUT_PATH / "growth_and_cue.csv")

    # ── Standardized PCA ──────────────────────────────────────────────────────
    print("\nStandardized PCA...")
    scaled_std, rxns_std, subs_std = prepare_matrix(flux_matrix)
    scores_std, loadings_std, ev_std = run_pca(scaled_std, subs_std, rxns_std)
    print(f"  Matrix {scaled_std.shape}  |  PC1+PC2 = {ev_std[0]+ev_std[1]:.1%} variance")

    scores_std.to_csv(OUT_PATH / "pca_scores_standardized.csv")
    loadings_std.to_csv(OUT_PATH / "pca_loadings_standardized.csv")

    plot_scores(
        scores_std, summary_df, ev_std,
        "MIT1002 flux PCA — standardized",
        OUT_PATH / "pca_scores_standardized.png",
    )
    plot_loadings(
        loadings_std, model, ev_std,
        "Standardized PCA",
        OUT_PATH / "loadings_standardized.png",
    )

    # ── Growth-rate-normalized PCA ────────────────────────────────────────────
    print("\nGrowth-rate-normalized PCA...")
    growth_rates = summary_df.loc[flux_matrix.index, "growth_rate"]
    scaled_gr, rxns_gr, subs_gr = prepare_matrix(flux_matrix, growth_rates=growth_rates)
    scores_gr, loadings_gr, ev_gr = run_pca(scaled_gr, subs_gr, rxns_gr)
    print(f"  Matrix {scaled_gr.shape}  |  PC1+PC2 = {ev_gr[0]+ev_gr[1]:.1%} variance")

    scores_gr.to_csv(OUT_PATH / "pca_scores_growth_normalized.csv")
    loadings_gr.to_csv(OUT_PATH / "pca_loadings_growth_normalized.csv")

    plot_scores(
        scores_gr, summary_df, ev_gr,
        "MIT1002 flux PCA — growth-rate-normalized",
        OUT_PATH / "pca_scores_growth_normalized.png",
    )
    plot_loadings(
        loadings_gr, model, ev_gr,
        "Growth-rate-normalized PCA",
        OUT_PATH / "loadings_growth_normalized.png",
    )

    # ── Growth and CUE summary ────────────────────────────────────────────────
    print("\nPlotting growth rate and CUE summary...")
    plot_growth_cue(summary_df, OUT_PATH / "growth_and_cue.png")

    print(f"\nAll results saved to: {OUT_PATH}")


if __name__ == "__main__":
    main()

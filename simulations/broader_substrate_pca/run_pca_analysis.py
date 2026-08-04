"""Broader substrate panel PCA of MIT1002 flux distributions.

Runs pFBA on a curated substrate panel at fixed total carbon uptake
(60 mmol C/gDW/hr), builds a (substrate x reaction) flux matrix, then runs
growth-rate-normalised PCA and saves scores coloured three ways:
  1. Chemical class (broad category)
  2. Metabolic entry point into central metabolism
  3. C:N ratio of the substrate (continuous, derived from model metabolite.elements)
"""

from pathlib import Path
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
from matplotlib.lines import Line2D
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

FILE_PATH = Path(__file__).resolve().parent
REPO_ROOT = FILE_PATH.parents[1]

import sys

sys.path.insert(0, str(REPO_ROOT))

from tools.media import MEDIA  # noqa: E402
DATA_DIR = REPO_ROOT / "data"
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)

TOTAL_UPTAKE = 60       # mmol C / gDW / hr
BIOMASS_RXN = "bio1_biomass"
CO2_EX_RXN = "EX_cpd00011_e0"

# ── Substrate annotations ─────────────────────────────────────────────────────

# Chemical class (broad category)
SUBSTRATE_CLASS = {
    "cpd00027": "Monosaccharide",   # Glucose
    "cpd00108": "Monosaccharide",   # Galactose
    "cpd00029": "Organic acid",     # Acetate
    "cpd00797": "Organic acid",     # 3-Hydroxybutyrate
    "cpd00123": "Organic acid",     # 3-Methyl-2-oxobutanoate (keto acid)
    "cpd00080": "Sugar phosphate",  # Glycerol-3-phosphate
    "cpd00035": "Amino acid",       # Alanine
    "cpd00023": "Amino acid",       # Glutamate
    "cpd00107": "Amino acid",       # Leucine
    "cpd00156": "Amino acid",       # Valine
    "cpd00322": "Amino acid",       # Isoleucine
    "cpd00051": "Amino acid",       # Arginine
    "cpd00039": "Amino acid",       # Lysine
    "cpd00069": "Amino acid",       # Tyrosine
    "cpd00129": "Amino acid",       # Proline
    "cpd00041": "Amino acid",       # Aspartate
    "cpd00033": "Amino acid",       # Glycine
    "cpd00127": "Other",            # Phenol
    "cpd23538": "Other",            # DHPS
}

CLASS_COLORS = {
    "Monosaccharide": "#2196F3",   # blue
    "Amino acid":     "#4CAF50",   # green
    "Organic acid":   "#F44336",   # red
    "Sugar phosphate": "#9C27B0",  # purple
    "Other":          "#9E9E9E",   # gray
}

# Metabolic entry point into central metabolism.
# Entry point = the central-metabolism intermediate the substrate first produces.
ENTRY_POINT = {
    # ── Glycolysis ──────────────────────────────────────────────────────────
    "cpd00027": "Glycolysis",          # Glucose      → G6P
    "cpd00108": "Glycolysis",          # Galactose    → G1P → G6P (Leloir pathway)
    "cpd00080": "Glycolysis",          # Glycerol-3-P → DHAP
    # ── Pyruvate ────────────────────────────────────────────────────────────
    "cpd00035": "Pyruvate",            # Alanine      → Pyr (alanine aminotransferase)
    "cpd00033": "Pyruvate",            # Glycine      → Ser → Pyr
    "cpd23538": "Pyruvate",            # DHPS         → C3 sulfo-intermediate → Pyr + sulfite
    # ── Acetyl-CoA ──────────────────────────────────────────────────────────
    "cpd00029": "Acetyl-CoA",          # Acetate      → AcCoA (acetyl-CoA synthetase)
    "cpd00797": "Acetyl-CoA",          # 3-HB         → AcAcCoA → 2× AcCoA
    "cpd00107": "Acetyl-CoA",          # Leucine      → HMG-CoA → AcCoA (ketogenic only)
    "cpd00039": "Acetyl-CoA",          # Lysine       → saccharopine/pipecolate → AcCoA
    # ── TCA — α-ketoglutarate ───────────────────────────────────────────────
    "cpd00023": "TCA — α-KG",          # Glutamate    → α-KG (glutamate dehydrogenase)
    "cpd00129": "TCA — α-KG",          # Proline      → Glu → α-KG
    "cpd00051": "TCA — α-KG",          # Arginine     → Glu → α-KG (succinyltransferase)
    # ── TCA — oxaloacetate ──────────────────────────────────────────────────
    "cpd00041": "TCA — OAA",           # Aspartate    → OAA (aspartate aminotransferase)
    # ── TCA — succinyl-CoA (via propionyl-CoA) ──────────────────────────────
    "cpd00156": "TCA — Succinyl-CoA",  # Valine       → propionyl-CoA → succinyl-CoA
    "cpd00322": "TCA — Succinyl-CoA",  # Isoleucine   → propionyl-CoA → succinyl-CoA (+AcCoA)
    "cpd00123": "TCA — Succinyl-CoA",  # KIC          → isobutyryl-CoA → propionyl-CoA → succinyl-CoA
    # ── Aromatic catabolism (split TCA entry) ───────────────────────────────
    "cpd00069": "Aromatic catabolism", # Tyrosine     → homogentisate → fumarate + AcAcCoA
    "cpd00127": "Aromatic catabolism", # Phenol       → catechol → β-ketoadipate → succinyl-CoA + AcCoA
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


# ── Helpers ───────────────────────────────────────────────────────────────────

def count_carbons(formula: str) -> Optional[int]:
    """Return number of carbon atoms from a molecular formula string."""
    if not formula:
        return None
    m = re.search(r"C(\d*)", formula)
    if m:
        n = m.group(1)
        return int(n) if n else 1
    return 0


def get_cn_ratio(model: cobra.Model, met_id: str) -> float:
    """Return C:N ratio for a substrate using model metabolite.elements.

    Tries the extracellular (_e0) then cytoplasmic (_c0) metabolite.
    Returns inf if nitrogen content is zero (no N in formula).
    Returns nan if the metabolite cannot be found at all.
    """
    for suffix in ("_e0", "_c0"):
        try:
            elems = model.metabolites.get_by_id(met_id + suffix).elements
            c = elems.get("C", 0)
            n = elems.get("N", 0)
            return float("inf") if n == 0 else c / n
        except KeyError:
            continue
    return float("nan")


# ── Substrate loading ─────────────────────────────────────────────────────────

def load_substrates(model: cobra.Model, media_defs: dict) -> pd.DataFrame:
    """Build substrate panel from the known-growth-phenotypes TSV.

    Filters for confirmed single-substrate growth, deduplicates by met_id,
    appends aspartate and glycine if absent, verifies exchange reactions exist.
    """
    tsv = DATA_DIR / "known_growth_phenotypes.tsv"
    df = pd.read_csv(tsv, sep="\t")

    df = df[df["growth"] == "Yes"].copy()
    df = df[~df["met_id"].astype(str).str.contains(",", na=True)].copy()

    # Normalise pro_exomet: fill missing as "No"
    if "pro_exomet" in df.columns:
        df["pro_exomet"] = df["pro_exomet"].fillna("No")
    else:
        df["pro_exomet"] = "No"

    # Sort so "Yes" pro_exomet rows come first — ensures dedup keeps the
    # annotated version when a substrate appears in multiple media contexts
    # (e.g. Valine: l1/NaN vs promm_no_c/Yes → keep promm_no_c row).
    df = (
        df.sort_values("pro_exomet", key=lambda s: s.map({"Yes": 0}).fillna(1),
                       ascending=True, kind="stable")
          .drop_duplicates(subset="met_id", keep="first")
          .copy()
    )

    # Manually add aspartate and glycine if absent from TSV.
    # Both were measured in the Pro exometabolome — mark as Yes.
    for name, met_id in [("Aspartate", "cpd00041"), ("Glycine", "cpd00033")]:
        if met_id not in df["met_id"].values:
            df = pd.concat(
                [df, pd.DataFrame([{"minimal_media": "l1", "c_source": name,
                                    "met_id": met_id, "growth": "Yes",
                                    "pro_exomet": "Yes"}])],
                ignore_index=True,
            )

    rxn_ids = {r.id for r in model.reactions}
    records = []

    for _, row in df.iterrows():
        met_id    = str(row["met_id"]).strip()
        c_source  = str(row["c_source"]).strip()
        ex_id     = f"EX_{met_id}_e0"
        media_key = str(row["minimal_media"]).strip()

        if ex_id not in rxn_ids:
            print(f"  SKIP (no exchange rxn)  : {c_source} ({met_id})")
            continue
        if media_key not in media_defs:
            print(f"  SKIP (unknown media '{media_key}'): {c_source}")
            continue

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
            "name":             c_source,
            "met_id":           met_id,
            "exchange_id":      ex_id,
            "media_key":        media_key,
            "n_c":              n_c,
            "substrate_class":  SUBSTRATE_CLASS.get(met_id, "Other"),
            "entry_point":      ENTRY_POINT.get(met_id, "Other"),
            "pro_exomet":       row.get("pro_exomet", "No"),
        })

    substrate_df = pd.DataFrame(records)
    print(f"\nSubstrate panel: {len(substrate_df)} substrates")
    print(
        substrate_df[["name", "met_id", "n_c", "media_key",
                       "substrate_class", "entry_point"]]
        .to_string(index=False)
    )
    return substrate_df


# ── pFBA simulations ──────────────────────────────────────────────────────────

def run_pfba(
    model: cobra.Model, media_defs: dict, substrate_df: pd.DataFrame
) -> tuple:
    """Run pFBA for every substrate; return (flux_matrix, summary_df)."""
    flux_records    = {}
    summary_records = []

    for _, row in substrate_df.iterrows():
        name      = row["name"]
        media     = media_utils.clean_media(model, media_defs[row["media_key"]])
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
                cue      = 1.0 - (co2_flux / TOTAL_UPTAKE)

                flux_records[name] = sol.fluxes.to_dict()
                summary_records.append({
                    "substrate":        name,
                    "met_id":           row["met_id"],
                    "growth_rate":      growth,
                    "co2_flux":         co2_flux,
                    "cue":              cue,
                    "substrate_class":  row["substrate_class"],
                    "entry_point":      row["entry_point"],
                    "pro_exomet":       row["pro_exomet"],
                })
                print(f"  OK   {name:28s}  mu={growth:.4f}  CUE={cue:.3f}")

            except Exception as exc:
                print(f"  FAIL {name}: {exc}")

    flux_matrix = pd.DataFrame(flux_records).T.fillna(0)
    summary_df  = pd.DataFrame(summary_records).set_index("substrate")
    return flux_matrix, summary_df


# ── PCA utilities ─────────────────────────────────────────────────────────────

def prepare_matrix(
    flux_matrix: pd.DataFrame,
    growth_rates: Optional[pd.Series] = None,
) -> tuple:
    """Normalise, drop zero/constant columns, standardise.

    Returns (scaled_array, reaction_names, substrate_names).
    """
    mat = flux_matrix.copy()
    if growth_rates is not None:
        mat = mat.div(growth_rates, axis=0)

    mat    = mat.loc[:, (mat != 0).any(axis=0)]
    scaled = StandardScaler().fit_transform(mat.values)

    valid     = ~np.isnan(scaled).any(axis=0)
    scaled    = scaled[:, valid]
    reactions = mat.columns[valid].tolist()
    return scaled, reactions, mat.index.tolist()


def run_pca(
    scaled: np.ndarray, substrates: list, reactions: list, n_components: int = 5
) -> tuple:
    """Fit PCA; return (scores_df, loadings_df, explained_variance_ratio)."""
    n   = min(n_components, scaled.shape[0] - 1, scaled.shape[1])
    pca = PCA(n_components=n)
    scores = pca.fit_transform(scaled)
    pcs    = [f"PC{i+1}" for i in range(n)]
    scores_df   = pd.DataFrame(scores, index=substrates, columns=pcs)
    loadings_df = pd.DataFrame(pca.components_.T, index=reactions, columns=pcs)
    return scores_df, loadings_df, pca.explained_variance_ratio_


# ── Plotting ──────────────────────────────────────────────────────────────────

def _annotate_points(ax, scores_df):
    """Add substrate name labels offset to the right of each point."""
    for sub, row in scores_df.iterrows():
        ax.annotate(
            sub, (row["PC1"], row["PC2"]),
            xytext=(7, 0), textcoords="offset points",
            fontsize=8, va="center", color="#444444",
        )


def _axis_labels(ax, ev):
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.axvline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel(f"PC1 ({ev[0]:.1%} variance)", fontsize=12)
    ax.set_ylabel(f"PC2 ({ev[1]:.1%} variance)", fontsize=12)


def plot_scores_categorical(
    scores_df: pd.DataFrame,
    color_series: pd.Series,
    color_dict: dict,
    legend_title: str,
    explained_var: np.ndarray,
    title: str,
    out_path: Path,
    pro_exomet_series: Optional[pd.Series] = None,
) -> None:
    """Scatter plot coloured by a categorical variable.

    If *pro_exomet_series* is provided, points with value "Yes" are drawn as
    stars (★) and all others as circles (●).  A second legend explains the
    marker shapes.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    for cls in sorted(color_series.unique()):
        color      = color_dict.get(cls, "#9E9E9E")
        class_mask = color_series == cls

        if pro_exomet_series is not None:
            # non-Pro → circle
            mask_other = class_mask & (pro_exomet_series != "Yes")
            subs_other = scores_df[mask_other]
            if len(subs_other):
                ax.scatter(subs_other["PC1"], subs_other["PC2"],
                           color=color, marker="o", s=180,
                           edgecolor="black", linewidth=0.6,
                           label="_nolegend_", zorder=3)
            # Pro exometabolite → star
            mask_pro = class_mask & (pro_exomet_series == "Yes")
            subs_pro = scores_df[mask_pro]
            if len(subs_pro):
                ax.scatter(subs_pro["PC1"], subs_pro["PC2"],
                           color=color, marker="*", s=300,
                           edgecolor="black", linewidth=0.6,
                           label="_nolegend_", zorder=3)
        else:
            subs = scores_df[class_mask]
            ax.scatter(subs["PC1"], subs["PC2"],
                       color=color, marker="o", s=180,
                       edgecolor="black", linewidth=0.6,
                       label="_nolegend_", zorder=3)

    _annotate_points(ax, scores_df)
    _axis_labels(ax, explained_var)
    ax.set_title(title, fontsize=13)

    # Legend 1: colour by class
    color_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=color_dict.get(cls, "#9E9E9E"),
               markeredgecolor="black", markersize=9, label=cls)
        for cls in sorted(color_series.unique())
    ]
    color_legend = ax.legend(handles=color_handles, title=legend_title,
                             bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.add_artist(color_legend)

    # Legend 2: marker shape (only when pro_exomet_series is provided)
    if pro_exomet_series is not None:
        shape_handles = [
            Line2D([0], [0], marker="*", color="w",
                   markerfacecolor="gray", markeredgecolor="black",
                   markersize=13, label="Pro exometabolite"),
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="gray", markeredgecolor="black",
                   markersize=9, label="Other"),
        ]
        ax.legend(handles=shape_handles, title="Source",
                  bbox_to_anchor=(1.02, 0), loc="lower left", fontsize=9)

    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_scores_continuous(
    scores_df: pd.DataFrame,
    value_series: pd.Series,
    cbar_label: str,
    explained_var: np.ndarray,
    title: str,
    out_path: Path,
    cmap: str = "plasma",
    inf_label: str = "No value (∞)",
    inf_color: str = "#aaaaaa",
    inf_marker: str = "^",
) -> None:
    """Scatter plot coloured by any continuous variable.

    Infinite/NaN values are plotted as grey triangles with a separate legend
    entry, so they don't compress the colourscale for finite values.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    finite_mask = np.isfinite(value_series)
    subs_fin = scores_df[finite_mask]
    vals     = value_series[finite_mask].values.astype(float)

    if len(vals) > 0:
        sc = ax.scatter(
            subs_fin["PC1"], subs_fin["PC2"],
            c=vals, cmap=cmap, vmin=vals.min(), vmax=vals.max(),
            s=180, edgecolor="black", linewidth=0.6, zorder=3,
        )
        cbar = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
        cbar.set_label(cbar_label, fontsize=10)

    subs_inf = scores_df[~finite_mask]
    if len(subs_inf) > 0:
        ax.scatter(
            subs_inf["PC1"], subs_inf["PC2"],
            color=inf_color, marker=inf_marker,
            s=180, edgecolor="black", linewidth=0.6,
            label=inf_label, zorder=3,
        )
        ax.legend(loc="lower right", fontsize=9)

    _annotate_points(ax, scores_df)
    _axis_labels(ax, explained_var)
    ax.set_title(title, fontsize=13)
    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


def plot_scores_cn_ratio(scores_df, cn_series, explained_var, title, out_path):
    """Convenience wrapper: colour by C:N ratio."""
    plot_scores_continuous(
        scores_df, cn_series,
        cbar_label="C:N ratio",
        explained_var=explained_var, title=title, out_path=out_path,
        cmap="plasma",
        inf_label="No nitrogen (C:N = ∞)", inf_color="#aaaaaa", inf_marker="^",
    )


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
        top_idx  = loadings_df[pc].abs().nlargest(top_n).index
        loadings = loadings_df.loc[top_idx, pc].sort_values()

        labels = []
        for rxn_id in loadings.index:
            try:
                rxn   = model.reactions.get_by_id(rxn_id)
                raw   = rxn.name or rxn_id
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


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading model...")
    model = cobra.io.read_sbml_model(REPO_ROOT / "model.xml")

    print("Loading media definitions...")
    media_defs = MEDIA

    print("\nBuilding substrate panel...")
    substrate_df = load_substrates(model, media_defs)

    print("\nRunning pFBA simulations...")
    flux_matrix, summary_df = run_pfba(model, media_defs, substrate_df)

    n_ok = len(summary_df)
    print(f"\nSuccessful: {n_ok}/{len(substrate_df)} substrates")
    if n_ok < 2:
        raise RuntimeError("Too few successful simulations to run PCA.")

    flux_matrix.to_csv(OUT_PATH / "flux_matrix.csv")
    summary_df.to_csv(OUT_PATH / "growth_and_cue.csv")

    # ── Active reaction counts ────────────────────────────────────────────────
    # Count reactions with |flux| > threshold per substrate.
    # Uses the raw (non-normalised) flux matrix so the count reflects total
    # metabolic activity, not per-unit-growth activity.
    ACTIVE_THRESHOLD = 1e-6
    active_rxn_counts = (flux_matrix.abs() > ACTIVE_THRESHOLD).sum(axis=1)
    active_rxn_counts.name = "active_reactions"
    print("\nActive reaction counts (|flux| > 1e-6):")
    for sub, n in active_rxn_counts.sort_values().items():
        print(f"  {sub:30s}  {n}")

    # ── C:N ratios from model metabolite.elements ─────────────────────────────
    cn_series = pd.Series(
        {row["substrate"]: get_cn_ratio(model, row["met_id"])
         for _, row in summary_df.reset_index().iterrows()},
        name="cn_ratio",
    )
    print("\nC:N ratios:")
    for sub, val in cn_series.items():
        label = f"{val:.2f}" if np.isfinite(val) else "∞ (no N)"
        print(f"  {sub:30s}  {label}")

    # ── Growth-rate-normalised PCA ────────────────────────────────────────────
    print("\nGrowth-rate-normalised PCA...")
    growth_rates = summary_df["growth_rate"]
    scaled, rxns, subs = prepare_matrix(flux_matrix, growth_rates=growth_rates)
    scores_df, loadings_df, ev = run_pca(scaled, subs, rxns)
    print(f"  Matrix {scaled.shape}  |  PC1+PC2 = {ev[0]+ev[1]:.1%} variance")

    scores_df.to_csv(OUT_PATH / "pca_scores_growth_normalized.csv")
    loadings_df.to_csv(OUT_PATH / "pca_loadings_growth_normalized.csv")

    # Colour by substrate class (★ = Pro exometabolite, ● = other)
    plot_scores_categorical(
        scores_df,
        color_series=summary_df.loc[scores_df.index, "substrate_class"],
        color_dict=CLASS_COLORS,
        legend_title="Chemical class",
        explained_var=ev,
        title="MIT1002 flux PCA (growth-rate-normalised) — by chemical class",
        out_path=OUT_PATH / "pca_scores_by_class.png",
        pro_exomet_series=summary_df.loc[scores_df.index, "pro_exomet"],
    )

    # Colour by entry point (★ = Pro exometabolite, ● = other)
    plot_scores_categorical(
        scores_df,
        color_series=summary_df.loc[scores_df.index, "entry_point"],
        color_dict=ENTRY_POINT_COLORS,
        legend_title="Entry point",
        explained_var=ev,
        title="MIT1002 flux PCA (growth-rate-normalised) — by entry point",
        out_path=OUT_PATH / "pca_scores_by_entry_point.png",
        pro_exomet_series=summary_df.loc[scores_df.index, "pro_exomet"],
    )

    # Colour by number of active reactions
    plot_scores_continuous(
        scores_df,
        value_series=active_rxn_counts.loc[scores_df.index],
        cbar_label="Active reactions (|flux| > 1e-6)",
        explained_var=ev,
        title="MIT1002 flux PCA (growth-rate-normalised) — by active reaction count",
        out_path=OUT_PATH / "pca_scores_by_active_reactions.png",
        cmap="viridis",
    )

    # Colour by C:N ratio
    plot_scores_cn_ratio(
        scores_df,
        cn_series=cn_series.loc[scores_df.index],
        explained_var=ev,
        title="MIT1002 flux PCA (growth-rate-normalised) — by C:N ratio",
        out_path=OUT_PATH / "pca_scores_by_cn_ratio.png",
    )

    plot_loadings(
        loadings_df, model, ev,
        "Growth-rate-normalised PCA",
        OUT_PATH / "loadings_growth_normalized.png",
    )

    print(f"\nAll results saved to: {OUT_PATH}")


if __name__ == "__main__":
    main()

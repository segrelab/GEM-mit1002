"""Run a PCA of full flux distributions across Pro top-10 substrate conditions
and save the results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
OUT_PATH.mkdir(exist_ok=True)


# Load the flux solutions stored in your single_and_cocktail_results.csv
# Each row is a condition; the 'fluxes' column is a dict of {rxn_id: flux}
results = pd.read_csv(
    FILE_PATH.parent
    / "single_and_cocktail_sims"
    / "results"
    / "single_and_cocktail_results.csv",
    converters={"fluxes": eval},
)

# Build a (conditions × reactions) matrix from the flux dicts
flux_matrix = pd.DataFrame(results["fluxes"].tolist(), index=results["substrate"])
flux_matrix = flux_matrix.fillna(0)

# Drop reactions that are zero in all conditions (no information)
flux_matrix = flux_matrix.loc[:, (flux_matrix != 0).any(axis=0)]
print(f"Matrix shape: {flux_matrix.shape}")

# Standardize: each reaction (column) gets mean=0, std=1
# This prevents high-flux reactions (e.g., ATP synthase) from dominating
scaler = StandardScaler()
flux_scaled = scaler.fit_transform(flux_matrix.values)

# Drop columns that became NaN after scaling (constant across conditions)
valid_cols = ~np.isnan(flux_scaled).any(axis=0)
flux_scaled = flux_scaled[:, valid_cols]
valid_reactions = flux_matrix.columns[valid_cols]
print(f"After dropping constants: {flux_scaled.shape}")

# Run PCA
pca = PCA(n_components=min(5, flux_scaled.shape[0] - 1))
scores = pca.fit_transform(flux_scaled)
explained_var = pca.explained_variance_ratio_

print(f"Variance explained per PC: {explained_var}")
print(f"Cumulative: {np.cumsum(explained_var)}")

# Save scores and loadings
scores_df = pd.DataFrame(
    scores,
    index=flux_matrix.index,
    columns=[f"PC{i+1}" for i in range(scores.shape[1])],
)
scores_df.to_csv(OUT_PATH / "pca_scores.csv")

loadings_df = pd.DataFrame(
    pca.components_.T,
    index=valid_reactions,
    columns=[f"PC{i+1}" for i in range(pca.n_components_)],
)
loadings_df.to_csv(OUT_PATH / "pca_loadings.csv")

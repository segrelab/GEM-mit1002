"""
Filter Katie's ProDiel exometabolite data to the 12 QA/QC'd metabolites
from Table 4.1 of her thesis chapter.

Inputs:
    ProDiel_quant_20260211.csv: full dataset with 84 metabolites,
        multiple analytical methods, intracellular + extracellular, 3 replicates,
        24 timepoints (0-46 h, every 2 h).

Outputs:
    ProDiel_filtered_replicates.csv: long-format, replicates preserved
    ProDiel_filtered_meanByTimepoint.csv: long-format, mean +/- SD across replicates
    ProDiel_filtered_wide.csv: wide-format (timepoint x metabolite), means only

Filtering logic:
    1. Keep only the 12 metabolites listed in Table 4.1 of the chapter
    2. Keep only extracellular measurements (INorEX == "EX")
    3. Set negative concentrations (from blank correction) to zero
    4. Average across triplicates at each timepoint
"""

from pathlib import Path

import pandas as pd

# --- Configuration ---

INPUT_FILE = Path("ProDiel_quant_20260211.csv")
OUTPUT_DIR = Path(".")

# Mapping from CSV CleanName -> chapter Table 4.1 name (for documentation)
TABLE_4_1_METABOLITES = {
    "2_keto_3_deoxy_6_phosphogluconate": "2-keto-3-deoxy-6-phosphogluconic acid",
    "3-phosphoglycerate": "3-phosphoglyceric acid",
    "adenosine_5_monophosphate": "5'-adenosine monophosphate",
    "uridine_5_monophosphate": "5'-uridine monophosphate",
    "a-ketoglutarate": "a-ketoglutaric acid",
    "cysteate": "cysteic acid",
    "desthiobiotin": "desthiobiotin",
    "glutamic_acid": "glutamic acid",
    "glutathione_oxidized": "glutathione (oxidized)",
    "glycerol_3_phosphate": "glycerol-3-phosphate",
    "n_acetylglutamic_acid": "N-acetyl-L-glutamic acid",
    "phosphoenolpyruvate": "phosphoenolpyruvic acid",
}


def filter_prodiel_data(input_file: Path, output_dir: Path) -> None:
    """Filter the full ProDiel dataset and write three output files."""
    df = pd.read_csv(input_file)

    # Confirm each Table 4.1 metabolite is present in the data
    csv_names = set(df["CleanName"].unique())
    print("Checking presence of Table 4.1 metabolites in input:")
    for csv_name, paper_name in TABLE_4_1_METABOLITES.items():
        present = csv_name in csv_names
        marker = "[OK]" if present else "[MISSING]"
        print(f"  {marker} {csv_name:40s}  ({paper_name})")

    # Filter to the target metabolites and extracellular measurements only
    df_filt = df[
        (df["CleanName"].isin(TABLE_4_1_METABOLITES.keys())) & (df["INorEX"] == "EX")
    ].copy()
    print(f"\nRows after metabolite + EX filter: {len(df_filt)}")
    print(f"Methods present: {df_filt['Method'].value_counts().to_dict()}")

    # I could only use metabolites that were measured by the aniline method
    # But that would cause me to lose 3-phosphoglycerate
    # To do that use the following
    # Restrict to aniline-derivatized method to match chapter's reported values
    # df_filt = df_filt[df_filt["Method"] == "AN"].copy()
    # print(f"Rows after restricting to Method=='AN': {len(df_filt)}")
    # dropped = (
    #     set(TABLE_4_1_METABOLITES) - set(df_filt["CleanName"].unique())
    # )
    # if dropped:
    #     print(
    #         f"WARNING: dropped by AN-only filter (not measured by aniline method): "
    #         f"{sorted(dropped)}"
    #     )

    # Handle negative concentrations from blank correction
    n_neg = (df_filt["nM"] < 0).sum()
    print(f"\nNegative nM values set to 0: {n_neg}")
    df_filt.loc[df_filt["nM"] < 0, "nM"] = 0

    # Aggregate across replicates
    agg = (
        df_filt.groupby(["CleanName", "timepoint"])
        .agg(
            mean_nM=("nM", "mean"),
            sd_nM=("nM", "std"),
            n_reps=("nM", "count"),
        )
        .reset_index()
    )
    print(f"\nFinal aggregated rows: {len(agg)}")
    print(f"Metabolites retained: {sorted(agg['CleanName'].unique())}")
    print(f"Timepoints: {sorted(agg['timepoint'].unique())}")

    # Write outputs
    output_dir.mkdir(parents=True, exist_ok=True)

    reps_path = output_dir / "ProDiel_filtered_replicates.csv"
    df_filt[["CleanName", "Method", "replicateID", "timepoint", "nM"]].to_csv(
        reps_path, index=False
    )
    print(f"\nWrote {reps_path}")

    mean_path = output_dir / "ProDiel_filtered_meanByTimepoint.csv"
    agg.to_csv(mean_path, index=False)
    print(f"Wrote {mean_path}")

    wide = agg.pivot(index="timepoint", columns="CleanName", values="mean_nM")
    wide_path = output_dir / "ProDiel_filtered_wide.csv"
    wide.to_csv(wide_path)
    print(f"Wrote {wide_path}")


if __name__ == "__main__":
    filter_prodiel_data(INPUT_FILE, OUTPUT_DIR)

"""
Plot results from the SMF dissipation sweep.

Key findings from run_smf_sweep.py:
  - Growth declines LINEARLY on all substrates (no threshold in growth curve)
  - The mechanism is PMF depletion: flagella short-circuits Na+ gradient, Na+/H+
    antiporter reverses to export excess Na+, consuming PMF, reducing ATP-synthase flux
  - NaNQR barely increases because the antiporter handles Na+ re-export instead
  - Acetate overflow INCREASES linearly with dissipation on most substrates
    (cell compensates for reduced ATP-synthase by using acetate kinase for SLP-ATP)
  - Substrate-specific onset: aspartate starts overflow at lb~76 (vs. baseline overflow
    on glucose/glycerol/alanine/lysine/glutamate)
  - Glycine: no overflow despite high Na+ symporter flux (routes carbons to CO2)

Figures produced:
  growth_vs_dissipation.png     -- growth rate across substrates
  overflow_vs_dissipation.png   -- acetate & NH3 secretion per substrate
  mechanism_vs_dissipation.png  -- ATP synthase and Na+/H+ antiporter fluxes (the why)
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams.update({"font.size": 9, "axes.linewidth": 0.8})

FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
CSV_PATH = OUT_PATH / "smf_sweep_results.csv"

# ---------------------------------------------------------------------------
# Load results
# ---------------------------------------------------------------------------
df = pd.read_csv(CSV_PATH)
feasible = df[df["feasible"]].copy()

# Substrate display order and styling
# Non-Na+ controls: dashed/dotted grays; Na+-symporter substrates: solid warm colors
SUBSTRATE_ORDER = ["glucose", "glycerol", "glutamate", "aspartate", "alanine", "glycine", "lysine"]

COLORS = {
    "glucose":   "#4d4d4d",
    "glycerol":  "#888888",
    "glutamate": "#d62728",
    "aspartate": "#ff7f0e",
    "alanine":   "#e377c2",
    "glycine":   "#9467bd",
    "lysine":    "#8c564b",
}
LINE_STYLES = {
    "glucose":   "--",
    "glycerol":  ":",
    "glutamate": "-",
    "aspartate": "-",
    "alanine":   "-",
    "glycine":   "-",
    "lysine":    "-",
}
NA_SYMPORTER = {
    "glucose":   False,
    "glycerol":  False,
    "glutamate": True,
    "aspartate": True,
    "alanine":   True,
    "glycine":   True,
    "lysine":    True,
}

# ---------------------------------------------------------------------------
# Figure 1: Growth vs. dissipation  (2 panels: absolute & normalised)
# ---------------------------------------------------------------------------
fig1, (ax_abs, ax_norm) = plt.subplots(1, 2, figsize=(9, 4))

for sub in SUBSTRATE_ORDER:
    sub_df = feasible[feasible["substrate"] == sub].sort_values("dissipation_lb")
    if sub_df.empty:
        continue
    g0 = sub_df["growth_rate"].iloc[0]
    kw = dict(color=COLORS[sub], ls=LINE_STYLES[sub], lw=1.6,
              label=sub, marker="o", ms=2.5, markevery=3)
    ax_abs.plot(sub_df["dissipation_lb"], sub_df["growth_rate"], **kw)
    ax_norm.plot(sub_df["dissipation_lb"], sub_df["growth_rate"] / g0, **kw)

for ax in (ax_abs, ax_norm):
    ax.set_xlabel("Forced Na⁺ import — SMF dissipation (mmol/gDW/hr)")
ax_abs.set_ylabel("Growth rate (hr⁻¹)")
ax_abs.set_title("Growth rate vs. SMF dissipation")
ax_norm.set_ylabel("Normalised growth (fraction of baseline)")
ax_norm.set_title("Growth rate vs. SMF dissipation (normalised)")
ax_norm.axhline(1.0, color="gray", lw=0.5, ls=":", zorder=0)

# Shared legend (right panel only)
ax_norm.legend(fontsize=7.5, loc="upper right",
               title="— Na⁺-symporter substrate\n-- glucose / glycerol (no Na⁺-symport)")

# Annotation: all slopes are similar (linear, substrate-independent)
ax_norm.text(
    0.05, 0.12, "All substrates: near-identical\nlinear decline\n(PMF depletion mechanism)",
    transform=ax_norm.transAxes, fontsize=7, va="bottom", color="gray",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="lightgray", alpha=0.8),
)

fig1.tight_layout()
fig1.savefig(OUT_PATH / "growth_vs_dissipation.png", dpi=150, bbox_inches="tight")
print("Saved: growth_vs_dissipation.png")

# ---------------------------------------------------------------------------
# Figure 2: Overflow — acetate + NH3 secretion per substrate
# ---------------------------------------------------------------------------
ncols = 4
nrows = 2
fig2, axes2 = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows), sharey=False)
axes2 = axes2.flatten()

for i, sub in enumerate(SUBSTRATE_ORDER):
    ax = axes2[i]
    sub_df = feasible[feasible["substrate"] == sub].sort_values("dissipation_lb")
    if sub_df.empty:
        ax.set_visible(False)
        continue

    x = sub_df["dissipation_lb"].values
    ace = sub_df["acetate"].fillna(0).values
    nh3 = sub_df["NH3_secretion"].fillna(0).values

    # Acetate (left y-axis)
    lace, = ax.plot(x, ace, color="#1f77b4", lw=1.5, label="acetate secretion")
    ax.fill_between(x, 0, ace, alpha=0.15, color="#1f77b4")
    ax.set_ylabel("Acetate secretion\n(mmol/gDW/hr)", color="#1f77b4", fontsize=7.5)
    ax.tick_params(axis="y", labelcolor="#1f77b4")

    # NH3 (right y-axis) — only non-negative values are meaningful secretion
    ax2 = ax.twinx()
    nh3_sec = np.clip(nh3, 0, None)   # uptake (negative) is not overflow
    lnh3, = ax2.plot(x, nh3_sec, color="#d62728", lw=1.5, ls="--", label="NH₃ secretion")
    ax2.set_ylabel("NH₃ secretion\n(mmol/gDW/hr)", color="#d62728", fontsize=7.5)
    ax2.tick_params(axis="y", labelcolor="#d62728")

    title_color = COLORS.get(sub, "black")
    na_label = " (Na⁺-symport)" if NA_SYMPORTER.get(sub) else " (no Na⁺-symport)"
    ax.set_title(f"{sub}{na_label}", color=title_color, fontweight="bold", fontsize=8.5)
    ax.set_xlabel("SMF dissipation (mmol/gDW/hr)", fontsize=7.5)

    # Annotate overflow onset for aspartate
    if sub == "aspartate" and (ace > 0.1).any():
        onset_lb = x[ace > 0.1][0]
        ax.axvline(onset_lb, color="gray", lw=0.8, ls=":", alpha=0.7)
        ax.text(onset_lb + 1, ax.get_ylim()[0], f"onset\nlb={onset_lb}",
                fontsize=6, color="gray", va="bottom")

# Hide the last unused subplot
axes2[-1].set_visible(False)

fig2.suptitle("Overflow secretion vs. SMF dissipation\n"
              "Acetate increases on all substrates; onset is substrate-specific",
              y=1.01, fontsize=9)
fig2.tight_layout()
fig2.savefig(OUT_PATH / "overflow_vs_dissipation.png", dpi=150, bbox_inches="tight")
print("Saved: overflow_vs_dissipation.png")

# ---------------------------------------------------------------------------
# Figure 3: Mechanism — ATP synthase and Na+/H+ antiporter reversal
# ---------------------------------------------------------------------------
fig3, (ax_atp, ax_anti) = plt.subplots(1, 2, figsize=(9, 4))

for sub in SUBSTRATE_ORDER:
    sub_df = feasible[feasible["substrate"] == sub].sort_values("dissipation_lb")
    if sub_df.empty:
        continue
    kw = dict(color=COLORS[sub], ls=LINE_STYLES[sub], lw=1.5,
              label=sub, marker="o", ms=2.5, markevery=3)
    ax_atp.plot(sub_df["dissipation_lb"], sub_df["ATPsyn_flux"], **kw)
    ax_anti.plot(sub_df["dissipation_lb"], sub_df["NaH_anti_flux"], **kw)

ax_atp.set_xlabel("SMF dissipation (mmol/gDW/hr)")
ax_atp.set_ylabel("ATP synthase flux (mmol/gDW/hr)\n(positive = ATP synthesis)")
ax_atp.set_title("ATP synthase flux\n(PMF depletion reduces ATP yield)")
ax_atp.legend(fontsize=7, loc="upper right")

ax_anti.axhline(0, color="gray", lw=0.8, ls=":", zorder=0)
ax_anti.set_xlabel("SMF dissipation (mmol/gDW/hr)")
ax_anti.set_ylabel("Na⁺/H⁺ antiporter flux (mmol/gDW/hr)\n(+ = Na⁺ in / H⁺ out   − = Na⁺ out / H⁺ in)")
ax_anti.set_title("Na⁺/H⁺ antiporter reversal\n(re-exports Na⁺ via PMF, not NADH)")
ax_anti.text(
    0.05, 0.12,
    "Antiporter reverses as SMF builds:\n"
    "Na⁺_c0 → Na⁺_e0 at cost of PMF.\n"
    "NaNQR (NADH cost) barely changes.",
    transform=ax_anti.transAxes, fontsize=7, va="bottom",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="lightgray", alpha=0.8),
)
ax_anti.legend(fontsize=7, loc="lower left")

fig3.tight_layout()
fig3.savefig(OUT_PATH / "mechanism_vs_dissipation.png", dpi=150, bbox_inches="tight")
print("Saved: mechanism_vs_dissipation.png")

# ---------------------------------------------------------------------------
# Figure 4: Combined summary — acetate overflow per substrate on one axis
# ---------------------------------------------------------------------------
fig4, ax4 = plt.subplots(figsize=(7, 4.5))

for sub in SUBSTRATE_ORDER:
    sub_df = feasible[feasible["substrate"] == sub].sort_values("dissipation_lb")
    if sub_df.empty:
        continue
    ace = sub_df["acetate"].fillna(0).values
    ax4.plot(
        sub_df["dissipation_lb"].values, ace,
        color=COLORS[sub], ls=LINE_STYLES[sub], lw=1.7,
        label=sub, marker="o", ms=3, markevery=3,
    )

ax4.set_xlabel("SMF dissipation (forced Na⁺ import flux, mmol/gDW/hr)")
ax4.set_ylabel("Acetate secretion (mmol/gDW/hr)")
ax4.set_title("Acetate overflow vs. SMF dissipation — all substrates")
ax4.axhline(0, color="black", lw=0.6)
ax4.legend(fontsize=8, loc="upper left",
           title="Na⁺-symporter subs: solid\nNo Na⁺-symport: dashed/dotted")

# Annotate the aspartate onset
aspartate_df = feasible[feasible["substrate"] == "aspartate"].sort_values("dissipation_lb")
ace_asp = aspartate_df["acetate"].fillna(0).values
x_asp = aspartate_df["dissipation_lb"].values
if (ace_asp > 0.1).any():
    onset = x_asp[ace_asp > 0.1][0]
    ax4.annotate(
        f"Aspartate onset\nlb = {onset}",
        xy=(onset, 0.1),
        xytext=(onset - 20, 3.5),
        arrowprops=dict(arrowstyle="->", color="#ff7f0e", lw=1),
        fontsize=7.5, color="#ff7f0e",
    )

fig4.tight_layout()
fig4.savefig(OUT_PATH / "acetate_overflow_summary.png", dpi=150, bbox_inches="tight")
print("Saved: acetate_overflow_summary.png")

# ---------------------------------------------------------------------------
# Print text summary
# ---------------------------------------------------------------------------
print("\n=== Summary: Acetate secretion at lb=0 (baseline) and lb=100 ===")
print(f"{'Substrate':<12} {'lb=0':>8} {'lb=100':>8} {'Increase':>8} {'Onset':>8}")
print("-" * 50)
for sub in SUBSTRATE_ORDER:
    sub_df = feasible[feasible["substrate"] == sub].sort_values("dissipation_lb")
    if sub_df.empty:
        continue
    ace = sub_df["acetate"].fillna(0).values
    x   = sub_df["dissipation_lb"].values
    ace_at_0   = ace[x == 0][0]  if (x == 0).any()   else float("nan")
    ace_at_100 = ace[x == 100][0] if (x == 100).any() else float("nan")
    increase = ace_at_100 - ace_at_0
    onset_mask = ace > 0.1
    onset = f"lb={x[onset_mask][0]}" if onset_mask.any() else "none"
    print(f"{sub:<12} {ace_at_0:>8.3f} {ace_at_100:>8.3f} {increase:>8.3f} {onset:>8}")

print("\n=== Summary: Growth slope (per 10 mmol/gDW/hr dissipation) ===")
for sub in SUBSTRATE_ORDER:
    sub_df = feasible[feasible["substrate"] == sub].sort_values("dissipation_lb")
    if len(sub_df) < 2:
        continue
    x = sub_df["dissipation_lb"].values
    g = sub_df["growth_rate"].values
    slope = np.polyfit(x, g, 1)[0] * 10
    print(f"  {sub:<12}: {slope:+.4f} hr⁻¹ per 10 mmol/gDW/hr SMF dissipation")

print("\n=== Mechanism summary: NaNQR vs. ATP synthase changes over lb 0→100 ===")
print(f"{'Substrate':<12} {'ΔNANQR':>8} {'ΔATPsyn':>10} {'ΔAntiport':>12}")
print("-" * 45)
for sub in SUBSTRATE_ORDER:
    sub_df = feasible[feasible["substrate"] == sub].sort_values("dissipation_lb")
    if sub_df.empty:
        continue
    x = sub_df["dissipation_lb"].values
    def delta(col):
        at0   = sub_df[col].values[x == 0][0]   if (x == 0).any()   else float("nan")
        at100 = sub_df[col].values[x == 100][0]  if (x == 100).any() else float("nan")
        return at100 - at0
    print(f"{sub:<12} {delta('NaNQR_flux'):>+8.2f} {delta('ATPsyn_flux'):>+10.2f} {delta('NaH_anti_flux'):>+12.2f}")

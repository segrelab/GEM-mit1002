# NADPH demand and the ED pathway — findings

## Question
Does rising NADPH (redox) demand pull glycolytic flux into the ED pathway in
the MIT1002 GEM? If so, redox demand could help explain why MIT1002 retains ED.

## Method
Added an O2-coupled, mass/charge-balanced NADPH "redox maintenance" drain
(`NADPH + H+ + 0.5 O2 -> NADP+ + H2O`), forced increasing flux through it, and
ran pFBA (max biomass) on glucose (EMP-preferring) and galacturonate (ED-obligate
control) at O2 = 20 and 1000. Markers: ED = `rxn01477` (eda), EMP = `rxn00558`
(PFK). See `run_nadph_titration.py` / `plot_nadph_titration.py`.

## Result: ED is never used as a NADPH source
- **Base titration (pFBA):** as forced NADPH demand rises on glucose, ED flux
  stays at 0. The demand is met entirely by NADP-isocitrate dehydrogenase
  (`rxn01387`, NADP-IDH). oxPPP, malic enzyme, transhydrogenase stay at 0.
  Same at both O2 levels. (Figure: `results/nadph_source_breakdown.png`.)
- **Loopless FBA:** under `loopless_solution`, IDH's contribution drops from
  4.48 to 0.67 — i.e. much of the pFBA "NADPH supply" was riding on
  thermodynamically-infeasible cycles. ED still = 0.
- **Knock out the main valves (IDH + malic enzyme + transhydrogenase):** ED
  still = 0. The model instead produces NADPH via reversible reductases run in
  reverse (e.g. proline:NADP+ oxidoreductase `rxn00931` at flux ~17; glutamate
  dehydrogenase `rxn00184`; choline:NADP+ oxidoreductase `rxn12191`).
- **Force NADPH through the glucose-6-P oxidation branch:** the model uses the
  oxidative PPP (`rxn01115`, 6PGDH), NOT ED. This is the biologically meaningful
  comparison: oxPPP yields 2 NADPH per glucose-6-P vs ED's 1, so ED has no redox
  advantage even against its own sibling branch.

## Why a fully clean test isn't tractable here
65 of 90 NADP-coupled reactions in the model are reversible. Many can act as
backdoor NADPH sources when run in reverse. Blanket or targeted blocking of these
(in their production direction) collapses growth to 0, because several — notably
glutamate dehydrogenase — are entangled with essential biosynthesis and nitrogen
assimilation. A rigorous redox analysis would require curating the directionality
of these ~64 reactions individually from thermodynamic data; doing so bluntly
makes the model non-viable.

## Bottom line for the paper
Within the GEM, ED confers **no stoichiometric payoff on glucose** — neither for
growth rate (growth *declines* with forced ED use; see `forced_ed_tradeoffs`) nor
for NADPH supply (IDH and oxPPP dominate; ED is never selected). The robustly
model-supported explanation for ED's presence is **substrate access** — uronic
acids that enter at KDPG obligately require it (see
`tradeoff_on_different_substrates`). This rules out both growth- and
redox-optimization rationales for ED on hexoses, which aligns with the goal of
challenging the "ED = thermodynamic/protein-cost optimum" view.

## Caveat worth a sentence / model-quality note
The large set of reversible NADP-coupled reactions is itself a model-quality
issue: it lets NADPH be generated almost anywhere, which would confound *any*
redox-based flux analysis. Worth flagging for future curation.

# Data

Data used to build, constrain and evaluate the MIT1002 model. This directory is
required by [standard-GEM](https://github.com/MetabolicAtlas/standard-GEM),
which asks that it carry a README describing how it is organised.

## Contents

| Path | What it is |
| --- | --- |
| `known_growth_phenotypes.tsv` | Experimentally observed growth or no-growth for a carbon source in a given medium, with a literature or lab reference. The basis of `test/test_growth.py` and the growth report. |
| `media_sources/` | Primary documents for the growth media — published recipes and lab protocols. Provenance for the definitions in `tools/media.py`. |
| `deprecated_identifiers/` | Reactions and metabolites removed from the model, and why. See the README there. |

## `known_growth_phenotypes.tsv`

One row per observed condition. Alongside the medium, the substrate and the
observed `growth` call, two columns control how the row is used.

### `growth`

`Yes`, `No`, or `Unsure`. `Unsure` is a statement about the observation itself —
growth was seen but was weak, slow, or not reproducible on transfer. Those rows
are simulated but not scored.

### `exclude_reason`

Empty means the row is scored. A value means the row is simulated and shown, but
left out of sensitivity and specificity, because comparing the model against it
would not tell you anything about the model.

`reason` is a closed vocabulary so the file stays queryable and the categories
stay meaningful; `test/test_phenotype_data.py` fails on any value not listed
here, so adding a category is a deliberate act that touches this README too.

| Value | Use when |
| --- | --- |
| `control_failed` | The experiment's own control did not behave as required, so the condition cannot be interpreted. |
| `id_uncertain` | The compound could not be confidently mapped to a model metabolite. |
| `conflicting_reports` | The same condition was scored differently by different sources and the disagreement is unresolved. |
| `not_representable` | A trusted result that flux balance analysis cannot reproduce in principle — regulation, inhibition, or a kinetic effect. |

The split that matters is between the first three and the last. The first three
are problems with the *observation*: a better experiment could fix them and the
row might come back into scoring. `not_representable` is a trusted observation
that a stoichiometric model cannot express even in principle. That is a
permanent property of the formalism rather than a data problem, and those rows
are worth reporting as a result rather than quietly dropping. Collapsing all
four into a single "ignore" flag would lose that distinction.

**A missing exchange reaction is not a reason to exclude a row.** If the model
has no way to take the compound up, that is a prediction of no growth and it is
scored as one — for most of these compounds the genome was searched and no
candidate transporter was found, which makes the absence a finding rather than a
gap. `not_representable` is for something else: a trusted observation that FBA
cannot express *even with* the right transporters in place, because the
mechanism is regulatory or kinetic. Reach for it when adding the missing
reaction would not resolve the disagreement.

### Excluding a row is not the same as deleting it

Excluded rows stay in the file. The experiment happened and the result is real;
what is in question is whether it can be compared against a model. Deleting the
row loses that history and invites someone to re-add the condition later without
knowing it was already looked at — the same reasoning behind
`deprecated_identifiers/`.

### Currently excluded

The nine `marine_broth_wo_yeast_and_peptone_no_n` conditions were run with a
no-nitrogen control that grew: the strain grows on glucose with no nitrogen
source added, so growth cannot be attributed to the nitrogen source supplied.

Six of the nine are marked `control_failed`. The three succinate rows are
deliberately **not** excluded, because the failure mode is asymmetric — it makes
extra nitrogen available, which can produce a spurious *positive* but cannot
produce a spurious *negative*. Those three conditions had nitrogen and still did
not grow, so they remain valid evidence that succinate is not a usable carbon
source, which is what `rxn05654_c0` was deprecated on.

Other rows worth considering for exclusion, left alone for now: `homarine` in
`mbm` (the notes flag an uncertain metabolite ID) and `proline`, which is scored
`Yes` in L1, `Unsure` in MBM and `No` in marine broth.

## Expected mismatches

`test/test_files/expected_phenotype_mismatches.tsv` records the mismatches that
are known and accepted at the current state of curation. `test_growth.py` fails
on a mismatch that is not listed, and equally on a listed mismatch that has
started passing — so an accidental fix is surfaced rather than silently rotting
the baseline. Regenerate it deliberately with
`python scripts/update_phenotype_baseline.py`, never automatically.

## Media: definitions live in code, not here

The media *definitions* are Python dictionaries in
[`tools/media.py`](../tools/media.py), not data files:

```python
from tools.media import MEDIA

model.medium = MEDIA["minimal_glucose"]
```

They are code rather than data because they are compositional — `minimal_glucose`
is `minimal` plus glucose, `promm_no_c` is `promm` minus its carbon sources — and
expressing that in a flat data file needs an inheritance mechanism that buys
nothing here.

What lives in this directory is the media *provenance*: the published recipes and
lab protocols in `media_sources/` that the definitions were transcribed from. If
you change a medium in `tools/media.py`, the source document here is what a
reviewer should be able to check it against.

### There used to be a pickle

Until recently the definitions were also serialised to
`test/test_files/media/media_definitions.pkl` and loaded with `pickle.load` from
22 places. That was removed because it was a cache of literal dicts — no
computation was being saved, so it bought nothing while costing three things: it
was a tracked binary so media changes could not be reviewed in a diff (and media
determine which growth phenotypes pass); it could go stale if someone edited the
definitions without regenerating it; and `pickle.load` is version-fragile and
executes arbitrary code on load.

The KBase and CarveMe tables previously derived from these definitions, along
with the scripts that generated them, were deleted rather than moved — they were
one-off exports for external tools and had gone stale. Regenerating anything
similar means writing it fresh against `tools.media.MEDIA`.

## A note on two key mismatches

`biomass/check_producibility.py` looks up media under the names `"mbm_media"` and
`"l1_media"`, but the registry keys — and the values in the `minimal_media`
column of `known_growth_phenotypes.tsv` — are `"mbm"` and `"l1"`. That lookup
could not have been matching before this refactor either; those names were never
in the pickle. Flagged rather than changed, since fixing it alters behaviour.

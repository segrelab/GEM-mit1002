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

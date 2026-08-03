# Deprecated identifiers

This directory is the **source of truth** for reactions and metabolites that were
once in the MIT1002 model and have since been removed.

Its purpose is to stop the same identifier being re-added, re-litigated, or
hunted for by someone who finds it referenced in an older figure, script, or
paper. Removing something from a GEM is a curation decision with as much
information content as adding it, and that decision should be recorded somewhere
more durable than a pull request thread.

## Files

| File | Contents |
| --- | --- |
| `deprecated_reactions.tsv` | Reactions removed from the model |
| `deprecated_metabolites.tsv` | Metabolites (SBML species) removed from the model |

## Schema

Both files share the same columns.

| Column | Required | Description |
| --- | --- | --- |
| `id` | yes | The identifier as it appeared in the model, **without** the SBML `R_`/`M_` prefix (e.g. `rxn00196_c0`, `cpd00225_c0`). This matches what COBRApy reports. |
| `name` | no | The `name` attribute the entity had when it was removed. Purely for human readability. |
| `reason` | yes | One value from the controlled vocabulary below. |
| `replaced_by` | no | Where this function lives in the model now, if anywhere. Semicolon-separated if more than one. See the section below — the name is slightly misleading. |
| `pr` | no | Pull request number that removed it, as `#123`. The PR remains the long-form record of the reasoning. **Leave this blank** — CI fills it in, see below. |
| `date` | no | `YYYY-MM-DD` the removal landed. |
| `notes` | no | At most one short line. Anything longer belongs in the PR. |

### Controlled vocabulary for `reason`

| Value | Use when |
| --- | --- |
| `duplicate` | Another reaction/metabolite in the model already did this. Set `replaced_by`. |
| `id_changed` | The entity survives under a different identifier. Set `replaced_by`. |
| `no_genomic_evidence` | No gene in the MIT1002 genome supports it. |
| `no_experimental_evidence` | Contradicted by, or unsupported by, phenotype data. |
| `mass_imbalance` | Unbalanced and not fixable without inventing chemistry. |
| `infeasible_cycle` | Participated in an ATP- or redox-generating cycle. |
| `dead_end` | Could never carry flux; blocked in all media. |
| `orphaned` | Only present because a reaction that used it was removed. Mostly metabolites. |
| `erroneous_annotation` | The underlying annotation was simply wrong. |
| `out_of_scope` | Real biology, but not part of what this model is meant to represent. |
| `unknown` | Backfilled from history and the reason could not be recovered. Do not use for new removals. |

`reason` is a closed vocabulary on purpose: it keeps the file queryable and
testable, and it keeps prose in the PR where it belongs. `test_deprecated.py`
fails on any value not in this list, so adding a category is a deliberate act
that touches this README too.

### What `replaced_by` actually means

The name is imperfect and worth reading carefully, because in the most common
case nothing was replaced.

`replaced_by` does not describe an event that happened to the model. It answers a
reader's question: **"I found this identifier somewhere — what should I look at
now?"** Three quite different situations produce a value, and only one of them is
a replacement in the ordinary sense:

- **`duplicate`** — the surviving reaction *was already there*. Nothing was
  substituted for anything. The model always represented this chemistry; it
  simply did so twice, and now does so once. `replaced_by` points at the
  survivor. Example: `rxn00154_c0` was removed as a duplicate of
  `rxn00011_c0 + rxn02342_c0 + rxn01801_c0`, all of which predated the removal.
- **`id_changed`** — a genuine rename. Same entity, new label. This is the only
  case where "replaced by" reads naturally.
- **One entity becomes several.** `rxn02200_c0` was removed for being a less
  accurate duplicate of `rxn02503_c0` + `rxn02201_c0`. That is not an
  equivalence; the function is now covered by a combination of reactions.

What is consistent across all three is that the column records **where the
biology lives now**, not what was done to the record. Read it as "use these
instead" rather than "this was swapped out for that."

#### Format

Semicolon-separated bare identifiers. No commas, no plus signs, no brackets:

```
rxn00011_c0; rxn02342_c0; rxn01871_c0; rxn01241_c0
```

Semicolons because a comma is a hazard in tabular data, and because Human-GEM
already uses semicolons for its own multi-value cells (`MNXR100067;MNXR100069`).
No brackets because this is a TSV cell, not a Python literal.

Reading and writing is forgiving, so you can paste a value straight out of a
commit message — `rxn00011_c0 + rxn02342_c0`, `[a, b]`, or plain whitespace
separation all parse, and `R_`/`M_` prefixes are stripped. Anything constructed
through the helper is normalised to the canonical form on write. But a
hand-edited file is checked strictly: `test_deprecated.py` fails on a cell
containing a comma, plus sign or bracket. Without that check the failure mode is
confusing, because a comma-separated cell parses as one long token and you get
told the target does not resolve rather than that the separator is wrong.

#### One reaction can be covered by several

The list carries no notion of "all of these together" versus "any of these".
`rxn00154_c0` was a single lumped reaction standing in for the entire pyruvate
dehydrogenase complex, and was removed because the model already contained the
four individual steps — so its `replaced_by` names all four, and the conjunction
is implied. Encoding that distinction in the identifier list would make the
column much harder to parse for very little benefit, so put the nuance in
`notes` instead:

| id | reason | replaced_by | notes |
| --- | --- | --- | --- |
| `rxn00154_c0` | `duplicate` | `rxn00011_c0; rxn02342_c0; rxn01871_c0; rxn01241_c0` | lumped reaction for the whole PDH complex; model already had the four individual steps |

That makes the empty-versus-filled distinction the one that matters most to
anyone using the model:

| `replaced_by` | Meaning |
| --- | --- |
| empty | The model no longer represents this at all. It is gone, deliberately. |
| filled | The model still represents this, under the listed identifier(s). Look there. |

A downstream user who finds a missing identifier in an old script needs exactly
that distinction, which is why `test_deprecated.py` requires the column for
`duplicate` and `id_changed` (where a claim about another identifier is implied
and useless without it) and checks that the targets actually resolve.

For the record, neither Human-GEM nor `standard-GEM` has a column like this.
Human-GEM's deprecated-identifier files are purely external-database
cross-references (`metBiGGID`, `metKEGGID`, `rxnRheaID`, and so on) and record
neither a reason nor a replacement; `standard-GEM` does not mention
deprecated identifiers at all. So this column is local to this repo, and the
name was kept only because renaming a column is churn — not because any standard
requires it.

## How to remove something from the model

Do not hand-edit these files and do not hand-edit `model.xml` to delete an
entity. Use the helper, which does both halves in one step:

```python
from tools.deprecate import deprecate_reactions

deprecate_reactions(
    ["rxn00196_c0"],
    reason="no_genomic_evidence",
    notes="no candidate gene in the MIT1002 genome",
)
```

### You do not need to know your PR number

Note the absence of a `pr` argument above. You remove things on your branch
*before* you open the pull request, so the number does not exist yet. Leave it
blank and CI fills it in: the `Custom-CI` workflow runs
`python -m tools.deprecate stamp-pr "$PR_NUMBER"` on every pull request and
commits the result, exactly as it already stamps the PR number into
`scripts/results/README.md`.

Stamping only fills **blank** cells. If you already know the relevant number —
usually because the removal closes an issue — pass it and CI will leave it
alone:

```python
deprecate_reactions(["rxn08703_c0"], reason="no_genomic_evidence", pr="#316")
```

The helper removes the reaction from `model.xml`, appends a row here, and
cascades to any metabolite or gene that the removal orphaned — logging the
orphaned metabolites with `reason="orphaned"` and an empty `replaced_by`, since
the model genuinely no longer represents them. This matters because
`test/test_sbml.py` already fails on isolated metabolites and genes, so a removal
that does not cascade breaks CI.

Run `--help` on the module for the metabolite equivalent and for a dry-run flag:

```bash
python -m tools.deprecate --help
```

## Why the model file also carries this information

A separate TSV has one real weakness: someone who downloads only `model.xml`
loses it. So `scripts/export_model.py` mirrors the identifier lists into the
SBML model's `<notes>` element, as:

```xml
<notes>
  <html xmlns="http://www.w3.org/1999/xhtml">
    <p>DEPRECATED_REACTIONS: rxn00196_c0; rxn01032_c0; ...</p>
    <p>DEPRECATED_METABOLITES: cpd00225_c0; ...</p>
    <p>DEPRECATED_INFO: full table with reasons at https://github.com/C-CoMP-STC/GEM-mit1002/tree/main/data/deprecated_identifiers</p>
  </html>
</notes>
```

`<notes>` rather than a custom `<annotation>` section, for a specific and
verified reason. A custom annotation in its own XML namespace *is* valid SBML —
libSBML accepts it and `cobra.io.validate_sbml_model` reports zero errors — but
COBRApy only parses SBO terms and RDF/MIRIAM CV terms out of annotations, so it
drops any foreign-namespace block on read and does not write it back. The
annotation would vanish the first time anyone round-tripped the model through
COBRApy, including our own export script, silently.

Model `<notes>` do survive. COBRApy parses `<p>KEY: VALUE</p>` pairs into
`model.notes` and writes them back out, the result is stable across repeated
round-trips, and it survives the JSON export too. The cost is that notes are a
flat string-to-string map: no tables, one line per key. Hence the split — the
TSV holds the structured detail, the notes carry the bare identifier list plus a
pointer.

The notes block is **generated**, never hand-written. `export_model.py`
regenerates it from these TSVs on every push to `main`, so if a contributor's
COBRApy script strips it, the next export puts it back. `test_deprecated.py`
checks the two agree.

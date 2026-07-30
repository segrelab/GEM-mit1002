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
| `replaced_by` | no | The identifier that now serves this function, if any. Semicolon-separated if more than one. Empty means nothing replaced it. |
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
orphaned metabolites with `reason="orphaned"` and `replaced_by` pointing back at
nothing. This matters because `test/test_sbml.py` already fails on isolated
metabolites and genes, so a removal that does not cascade breaks CI.

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

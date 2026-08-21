# Curation process analysis

Code and data for figure 2 of the manuscript: how the model's agreement with the
known growth phenotypes moved over the course of curation.

## Files

| File | What it is |
| --- | --- |
| `run_tests_on_prs.py` | Fetches `model.xml` as of each merged PR and scores it. Writes the two CSVs below. |
| `phenotype_confusion_over_time.csv` | **Source of truth.** One row per PR: the confusion matrix, every unscored category, sensitivity, specificity, model size. |
| `growth_match_summary.csv` | Narrow view of the same numbers, read by the plotting script. |
| `plot_match_over_time.py` | Draws figure 2B from `growth_match_summary.csv`. |
| `match_over_time.png` | The rendered panel. |

The confusion matrices shown in figure 2A are two rows of
`phenotype_confusion_over_time.csv` — the first PR and the current release.
Read them off the file; do not count them by hand. Counting them by hand is what
produced the version that did not add up.

## Regenerating

```bash
python curation_process/run_tests_on_prs.py   # needs `gh` authenticated
python curation_process/plot_match_over_time.py
```

Both are safe to re-run. `run_tests_on_prs.py` is incremental: it only evaluates
PRs it has no result for, PRs whose stored result is `ERROR`, and PRs whose
stored result was scored under an older `SCORING_VERSION`. Set
`FORCE_RERUN = True` to re-evaluate everything regardless.

**Re-run it after editing `data/known_growth_phenotypes.tsv`.** The stored rows
were scored against the table as it stood when they were computed, so adding a
condition leaves the early points on the line answering a different question
from the late ones. `test/test_curation_timeline.py` fails when that has
happened, by comparing each row's `Conditions` count against the current file.

## What counts as a match

Scoring lives in `tools/phenotypes.py`, shared with `test/test_growth.py`, so
the figure and the CI test cannot disagree.

A **match** is a true positive or a true negative. Only two kinds of condition
are held out:

- `Unsure` experimental results — there is nothing to agree or disagree with.
- Rows with an `exclude_reason` — see `data/README.md`.

Everything else is scored, including conditions whose compound has no exchange
reaction in the model. For most of those the absence of a transporter is a
finding rather than a gap: the genome was searched, no candidate was found, and
the model reflects that. "No genomic evidence for uptake, therefore no growth"
is a real mechanistic prediction, and it is most of what a stoichiometric model
has to say about *failure* to grow.

It is scored on the model's actual solution, not by treating a missing exchange
as an automatic "No" — those are different, because a multi-compound condition
can still grow on the compounds that are present. `Methionine, Pyruvate` and
`Cystine, Pyruvate` are both in that state and the shortcut gets both of them
backwards. See the "Missing exchange reactions" section of
`tools/phenotypes.py`.

`No Uptake Route` in the CSV counts how many of the model's no-growth
predictions rest on an absent transporter. It is a subset of `Scored`, not
another bucket, so it does not enter the `Scored + Unsure + Excluded + Invalid
Solve = Conditions` sum. Quote it alongside specificity: a specificity carried
mostly by absent transporters is a different claim from one carried by network
structure, and a reader is entitled to know which it is.

The denominator plotted is `Interpretable`: conditions with a definite Yes/No
and no exclusion reason. It is a property of the phenotype table, not of the
model, so it is identical for every point on the line, which is what makes the
two ends comparable.

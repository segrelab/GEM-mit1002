[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://hgscott.github.io/mit1002-model)

# MIT1002-model

This repo contains the *Alteromonas macleodii* MIT1002 model, and code associated with its creation, curation, and testing.

The model was generated using GenBank genome for MIT1002 (Accession Number: NZ_JXRW01000001), accessible via KBase.
The narrative for generating the draft model, is available here: https://narrative.kbase.us/narrative/208605

This repo uses GtiHub actions to automatically test the model.
Upon every push, pull request, manual trigger:
1. A new MEMOTE report is generated, and saved as "index.html"
2. Run custom tests
    * Validate the SBML file
    * Test for growth on no carbon sources
    * Test known growth phenotypes, and regenerate the experimental vs predicted growth heatmap figure
    * Run the MEMOTE test to search for ATP generating cycles
    * Check that nothing in the deprecated identifier lists is back in the model
3. The model is exported to JSON and excel formats

Note: MACAW is **not** run as part of the action due to the longer run time of the dilution test.
To run MACAW use:
```
python run_macaw.py
```

## Repository layout

Three directories hold code, and the distinction between them is about *what the
code does*, not what it is about. Please put new code in the matching one.

| Directory | Contains | How it runs |
| --- | --- | --- |
| `test/` | Checks that assert something about the model and pass or fail | Automatically, via `pytest` in CI. A failure blocks the PR |
| `scripts/` | Code that generates an artifact for a person to look at — a table, a plot, an exported file. No pass/fail | Automatically in CI, writing to `scripts/results/` |
| `tools/` | Importable functions, and command-line utilities a curator runs deliberately to *change* the model | By hand, or imported by the above |

`tools/deprecate.py` is the current example of the third kind: you invoke it
yourself when you remove something, and `scripts/export_model.py` and
`test/test_deprecated.py` both import from it.

Two other directories hold code that is neither of these: `curation_process/`
analyses the history of the curation effort itself across past PRs, and
`biomass/`, `genome/`, `escher/` and similar hold the exploratory work behind
particular parts of the model.

Note that [standard-GEM](https://github.com/MetabolicAtlas/standard-GEM), which
yeast-GEM and Human-GEM follow, asks for a single `code/` directory instead. The
split above is a deliberate refinement of that; `code/` is also a poor Python
package name because it shadows a standard-library module.

## To contribute to the model
1. Make a GitHub account
2. Make a fork/branch of this repo
3. Make your edits to the model on the XML file
4. If you are *removing* a reaction or metabolite, use the deprecation helper rather than deleting it by hand (see below)
5. Open a pull request

## Removing reactions and metabolites

Reactions and metabolites that have been removed from the model are recorded in
[`data/deprecated_identifiers/`](data/deprecated_identifiers/). Removal is a
curation decision with as much information content as an addition, and recording
it stops the same identifier being re-added or hunted for by someone who found it
in an older figure or script.

Remove things with the helper, which edits the model and updates the list in one
step, and cleans up any metabolite or gene the removal orphaned:

```
python -m tools.deprecate reaction rxn00196_c0 \
    --reason no_genomic_evidence --dry-run
```

Drop `--dry-run` to actually apply it. `--reason` takes a fixed vocabulary
documented in
[`data/deprecated_identifiers/README.md`](data/deprecated_identifiers/README.md);
the full reasoning still belongs in the pull request description, which the list
links back to.

You do not need to pass a PR number — you do not have one yet when you are
working on your branch. CI fills it in on every pull request and commits the
result, the same way it stamps the PR number into `scripts/results/README.md`.

The identifier lists are also mirrored into the SBML model's `<notes>`, so a
person who downloads only `model.xml` can still tell that those identifiers were
deliberately removed and where to find the reasons. `test/test_deprecated.py`
fails if the model and the lists disagree.

## Setting Up the Environment
To ensure a smooth setup and avoid system conflicts, follow these steps to create and activate a Python virtual environment before installing dependencies.

1. Check Your Python Version
First, make sure you have Python 3.11 or 3.10 installed.
Run the following command to check your version:
```
python3 --version
```
If it shows Python 3.13, we strongly recommend using Python 3.11 or 3.10, as some dependencies may not yet support Python 3.13.

To install Python 3.11 via Homebrew (if needed), run:
```
brew install python@3.11
```
2. Create a Virtual Environment

Once you have the correct Python version, create a virtual environment:
```
python3.11 -m venv .venv  # Use python3.10 if needed
```
This creates a .venv folder in your project directory, isolating dependencies from the system Python.

3. Activate the Virtual Environment

Before installing packages, activate the virtual environment:

Mac/Linux:
```
source .venv/bin/activate
```
Windows (Command Prompt):
```
.venv\Scripts\activate
```
Windows (PowerShell):
```
.\.venv\Scripts\Activate
```
Your terminal should now show (.venv) at the beginning of the prompt, indicating the environment is active.

4. Upgrade Pip

Inside the virtual environment, upgrade pip to avoid compatibility issues:
```
pip install --upgrade pip setuptools
```
5. Install Dependencies

Now, install all required dependencies:
```
pip install -r requirements.txt
```
If you encounter an error about "externally managed environment", add the following flag:
```
pip install -r requirements.txt --break-system-packages
```
6. Verify Installation

To ensure everything is working correctly, run:
```
python --version  # Should be 3.11 or 3.10
pip list  # Should show installed dependencies
```
If everything looks good, you're ready to start using the project! 🎉

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

## To contribute to the model
1. Make a GitHub account
2. Make a fork/branch of this repo
3. Make your edits to the model on the XML file
4. If you are *removing* a reaction or metabolite, use the deprecation helper rather than deleting it by hand (see below)
5. Open a pull request

## Removing reactions and metabolites

Reactions and metabolites that have been removed from the model are recorded in
[`data/deprecatedIdentifiers/`](data/deprecatedIdentifiers/). Removal is a
curation decision with as much information content as an addition, and recording
it stops the same identifier being re-added or hunted for by someone who found it
in an older figure or script.

Remove things with the helper, which edits the model and updates the list in one
step, and cleans up any metabolite or gene the removal orphaned:

```
python -m scripts.deprecate reaction rxn00196_c0 \
    --reason no_genomic_evidence --pr '#317' --dry-run
```

Drop `--dry-run` to actually apply it. `--reason` takes a fixed vocabulary
documented in
[`data/deprecatedIdentifiers/README.md`](data/deprecatedIdentifiers/README.md);
the full reasoning still belongs in the pull request description, which the list
links back to.

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

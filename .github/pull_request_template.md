#### Main improvements in this PR:
<!-- Try to be as clear as possible in providing a summary of the work and reference the corresponding issue.
Is it fixing/adding something in the model?
Is it an additional test/function/dataset? 
e.g. This PR improves/fixes # by ...
-->


**I hereby confirm that I have:**
<!-- *Note: replace [ ] with [X] to check the box. -->
- [ ] Made my edits to the model on the XML file
- [ ] Tested my code on my own computer for running the model
- [ ] Selected `dev` as a target branch
- [ ] Removed any reactions/metabolites using `python -m scripts.deprecate`, so they are recorded in `data/deprecatedIdentifiers/` (see the [README](../data/deprecatedIdentifiers/README.md) there). Nothing to remove in this PR? Check the box.

<!-- The deprecate helper removes the entity and logs it in one step, and cleans
up metabolites/genes the removal orphaned. Doing it by hand is how the list goes
stale, and test_deprecated.py will fail if the model and the list disagree.

    python -m scripts.deprecate reaction rxn00196_c0 \
        --reason no_genomic_evidence --pr '#317' --dry-run

Use the PR description above for the full reasoning -- the TSV only records a
one-word reason category and a link back here. -->


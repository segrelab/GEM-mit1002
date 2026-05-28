# ProDiel Data

## ProDiel_quant_20260211.csv
This table includes both intra- and extracellular metabolite concentrations (for variable INorEX, IN = intracellular, EX = extracellular). We have a few stray QA/QC strands to wrap up (eg, some metabolites are measured through multiple methods – the dataset I’ve attached chooses just one method per metabolite, but we’ve got a little double checking to do).

## extrapolated_cellcounts.csv
This table includes the cell abundance data for Prochlorococcus over the diel cycle, measured by flow cytometry. The columns are:
* "comment": Seems to be a combination of code for the replicate (a1, b1, or c1) and the time point (e.g., t4, t32, etc.)
* "time_h": the hour (every other hour 2-48) of the time point
* "typeOfData": Was this data lost and extrapolated (for t<12), measured (for t >= 12), or interperolated (for t=26 and t=34, where one replicate was missing. That replicate was filled in to have the mean value of the other two.)
* "cell_count_mean": the mean across the three replicates at that timepoint. This is the same value repeated across the three rows at any given timepoint.
* "cell_count_sd": the standard deviation across the three replicates at that timepoint. Also repeated across the three rows at each timepoint.
* "cell count": the individual replicate's measured (or extrapolated/interpolated) cell density at that timepoint, in cells/mL
Samples before t=12 were lost due to an instrument error, and so were extrapolated. Katie did the extrapolation, so the data from this table will match her analyses.
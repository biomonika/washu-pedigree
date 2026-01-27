# Scripts for CDR inheritance tracing

In order to trace CDRs across the generations, three- and two- generational centromere assigments are taken from the file `verkko_2.0_centromere_transmissions_chromosome_transmissions.csv`. To reduce the effect of noise and outliers (occassionaly, small fluctuations/dips in methylation can generate a subCDR that is very distant from the main cluster), subCDRs smaller than 3kb were filtered out, and outliers were removed by the script based on the DBSCAN algorithm (`keep_main_dbscan_cluster_filter_outliers.py`). 

Finally, the comparisons between generations 1 and 2, as well as generations 2 and 3, are calculated using the script `compare_gens_overlap_span.R`. The inter-generational analysis includes subCDR overlap across all three generations, midpoint absolute distances acros generations, as well as inter-generationl differences in total subCDR (sum of all subCDRs) and CDR length (representing total span from the first basepair of the first subCDR to the last basepair of the last subCDR).

The whole workflow is orchestrated by the script `parse_transmissions.py`. 

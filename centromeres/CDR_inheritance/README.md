# Scripts for CDR inheritance and positioning analysis

The CDR positioning changes across generations were calculated by re-mapping the sequencing reads from all the individuals on the mother (PAN027) as a reference, allowing for a shared-coordinate system. The mapping was performed using minimap2 with the map-ont setting. 

## CDR inheritance across the 3 generations
In order to trace CDRs across the generations, three- and two- generational centromere assigments are taken from the file `verkko_2.0_centromere_transmissions_chromosome_transmissions.csv`. To reduce the effect of noise and outliers (occassionaly, small fluctuations/dips in methylation can generate a subCDR that is very distant from the main cluster), subCDRs smaller than 3kb were filtered out, and outliers were removed by the script based on the DBSCAN algorithm (`keep_main_dbscan_cluster_filter_outliers.py`). 

Finally, the comparisons between generations 1 and 2, as well as generations 2 and 3, are calculated using the script `compare_gens_overlap_span.R`. The inter-generational analysis includes subCDR overlap across all three generations, midpoint absolute distances acros generations, as well as inter-generationl differences in total subCDR (sum of all subCDRs) and CDR length (representing total span from the first basepair of the first subCDR to the last basepair of the last subCDR).

The whole workflow is orchestrated by the script `parse_transmissions.py`. 

## CDR inheritance across the 2 generations

In order to study the inheritance across two generations only (centromeres that were transmitted from grandparents to the mother bus subsequently not transmitted to the third generation), we used the main coordinating script `parse_transmissions_2G.py`. 

This script uses the transmission assignments from the file `2G.centromere.transmissions.txt` and calculates the variability in the positioning using the script `compare_gens_overlap_span_2G.R`. 

## CDR positioning concordance or discordance

CDRs can be in concordant or discordant positions within diploid genomes. This analysis is performed using the script `centromere_positioning.Rnw`


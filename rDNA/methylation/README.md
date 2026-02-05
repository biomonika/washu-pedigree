# pipeline for methylation calling
`arr_job.sh` calls `pipeline.sh` for all samples, pipeline calls everything else

First `map_to_ref_and_split_per_group.sh` maps rDNA ref to reads and handles file structure setup. This calls `filter_both.py` and `remove_chimeras.py` to filter low quality alignments and chimeric ONT reads.

`pipeline.sh` then calls `get_18S_cov.sh` which uses `filter_both_18S.py` to get coverage of the 18S gene for later copy num estimation

`pipeline.sh` then calls `get_methylation_whole_group.sh` to get methylation of rDNA unit for entire sample (or just of 45S gene if passing a subregion BED file) using `modkit`

`pipeline.sh` then calls `get_methylation_breakdown.sh` to get methylation of individual rDNA units on a per-read per-unit basis (or just of 45S gene if passing a subregion BED file) using `modkit`

`pipeline.sh` can then optionally call `get_TR_breakdown.sh`, but that is not relevant to the scope of this paper and the script is not included here

Then to summarize output:
`summary_and_analysis.sh` calls `get_read_summary.sh` which uses the above data to estimate methylation of individual rDNA units or subregion, and copy number of the 18S gene, as well as outputs a distribution of unit methylation using `get_sample_distribution.py`. `summary_and_analysis.sh` also optionally calls other scripts outside the scope of this paper and not included here.

`summarize_group_meth.sh` will calculate avg methylation over the region of interest for the entire group.

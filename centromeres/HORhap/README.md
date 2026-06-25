# HORhap analysis

## Generate HORhap annotations

HORhap annotations were generated using [horhap_tool](https://github.com/fedorrik/horhap_tool).

1. For each chromosome separately, HOR sequences were extracted from all five assemblies using the StV annotations and aligned with MUSCLE3.
2. A pairwise distance matrix of HORs was calculated, and Ward clustering was performed.
3. The clustering dendrogram was cut at a selected depth to define nine clades, hereafter referred to as HORhaps.
4. HORhap assignments were added to HORs in the StV annotation and colored according to HORhap.

## Align transmitted HORhap annotations

HORhap annotations of transmitted centromeres were aligned, and indel events were called using [annotaligner](https://github.com/fedorrik/annotaligner).

1. For each pair of transmitted centromeres, HORhap annotations were aligned using the Needleman-Wunsch algorithm, with the combined StV and HORhap assignment of each HOR used as the alignment token instead of nucleotides.
2. Gaps in the resulting alignments were interpreted as indel events.

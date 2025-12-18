## 1. Align pedigree graph nodes to assemblies
```
/minimap2 -x map-hifi -t24 assembly.v1.0.PAN027.diploid.fa assembly_threegenerational_nodes_maternal.fasta -o threegen_maternal_to_PAN027.v1.0.diploid.paf
```
```
minimap2 -x map-hifi -t24 assembly.v1.0.PAN027.diploid.fa assembly_threegenerational_nodes_paternal.fasta -o threegen_paternal_to_PAN027-assembly.v1.0.diploid.paf
```
```
minimap2 -x map-hifi -t24 assembly.v1.0.PAN027.diploid.fa assembly_twogenerational_nodes_maternal.fasta -o twogen_maternal_to_PAN027-assembly.v1.0.diploid.paf
```
```
minimap2 -x map-hifi -t24 assembly.v1.0.PAN027.diploid.fa assembly_twogenerational_nodes_paternal.fasta -o twogen_paternal_to_PAN027-assembly.v1.0.diploid.paf
```
## 2. Extract bed files from paf

```
python extract_regions_from_paf.py threegenerational_nodes_maternal_contigs.txt threegen_maternal_to_PAN027-assembly.v1.0.diploid.paf threegen_maternal_contigs_mappedregions_assembly.v1.0.diploid.bed
```
```
python extract_regions_from_paf.py threegenerational_nodes_paternal_contigs.txt threegen_paternal_to_PAN027-assembly.v1.0.diploid.paf threegen_maternal_contigs_mappedregions_assembly.v1.0.diploid.bed
```
## 3. Get CHM13 coordinates using RagTag scaffolding

```
ragtag.py scaffold chm13v2.0.fa.gz assembly.v1.0.PAN027.maternal.fa -o ragtag-maternal -t 24 --aligner minimap2 --mm2-params '-x map-hifi'
```
```
ragtag.py scaffold chm13v2.0.fa.gz assembly.v1.0.PAN027.paternal.fa -o ragtag-paternal -t 24 --aligner minimap2 --mm2-params '-x map-hifi'
```




# Segmental duplications 

Segmental duplications were annotated using biser v1.4 `run_biser.sh` with the settings `biser -o biser.${assembly_name} --keep-contigs -t 2 --gc-heap 10G ${assembly}`. The output in the bedpe format was additionally also formatted into a bed format (`process_biser.sh`). Biser requires softmasked assemblies as an input; therefore, each haplotype was first softmasked using repeatmasker annotations and then processed individually. 

The repeatmasker annotations were obtained as part of the cenSat worflow `https://github.com/kmiga/alphaAnnotation.git`, using commit `3a0bdd6`. 
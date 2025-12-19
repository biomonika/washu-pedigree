#!/bin/bash

# Lift over local ancestry results from CHM13 to sample-specific assembly coordinates

#================#
# Set data paths #
#================#

# Path to UCSC liftover binary
liftOver=/home/syan11/code/liftOver
# Path to directory containing chain files between CHM13 and sample assemblies
chainPath=/scratch16/rmccoy22/syan11/washu_local_ancestry/chain_files_pedigree_assemblies_v1.0
# Path to directory containing sample-specific fasta files
assemblyPath=/scratch16/rmccoy22/syan11/washu_local_ancestry/08_ancestry_plots/pedigree_assemblies_v1.0_censats
# Path to flare output local ancestry VCFs
flarePath=/scratch16/rmccoy22/syan11/washu_local_ancestry/flare
# Output directory
flareOutPath=flare_lifted

ml gatk


#===========================================================#
# Edit all chain files to remove assembly-specific prefixes #
#===========================================================#

for sample in PAN010 PAN011 PAN028
do
    for hap in haplotype1 haplotype2
    do
        chain=$chainPath/CHM13_to_assembly.v1.0.$sample.$hap.chain
        chainNew=$chainPath/CHM13_to_assembly.v1.0.$sample.$hap.edited.chain
        cat $chain | sed "s/$sample.//g" | sed "s/.$hap//g" > $chainNew
    done
done
# Do PAN027 separately since it has different haplotype names
for hap in maternal paternal
do
    chain=$chainPath/CHM13_to_assembly.v1.0.PAN027.$hap.chain
    chainNew=$chainPath/CHM13_to_assembly.v1.0.PAN027.$hap.edited.chain
    cat $chain | sed "s/PAN027.//g" | sed "s/.$hap//g" > $chainNew
done


#============================================================#
# Edit reference fastas to remove assembly-specific prefixes #
#============================================================#

for sample in PAN010 PAN011 PAN027 PAN028
do
    # Get names of assembly haplotypes
    hap1=haplotype1
    hap2=haplotype2
    if [ "$sample" == "PAN027" ]; then
        hap1=maternal
        hap2=paternal
    fi

    # Edit fasta files
    inFasta=`ls $assemblyPath/assembly.v1.0.${sample}.diploid-*.fa`
    echo $inFasta
    for hapName in $hap1 $hap2
    do
        outFasta=$assemblyPath/assembly.v1.0.${sample}.${hapName}.fa
        # Separate haplotypes and edit contig names
        grep -A1 $hapName $inFasta | \
            sed "s/$sample.//g" | sed "s/.$hapName//g" \
            > $outFasta
        # Create .dict file
        gatk CreateSequenceDictionary --REFERENCE $outFasta
    done
done


#===================================================#
# Liftover ancestry from CHM13 to assembly-specific #
#===================================================#

for sample in PAN011 PAN010 PAN027 PAN028
do
    # Get names of assembly haplotypes
    hap1=haplotype1
    hap2=haplotype2
    if [ "$sample" == "PAN027" ]; then
        hap1=maternal
        hap2=paternal
    fi

    for hapName in $hap1 $hap2
    do
        # Get reference fasta and chain file
        referenceFasta=$assemblyPath/assembly.v1.0.${sample}.${hapName}.fa
        chain=$chainPath/CHM13_to_assembly.v1.0.$sample.$hapName.edited.chain

        # Liftover
        for i in {1..22} X_nonpar X_par1 X_par2
        do
            flareVcf=$flarePath/${sample}_chm13_chr$i.anc.vcf.gz
            liftedVcf=$flareOutPath/$sample.chr$i.anc.${hapName}_lifted.vcf.gz
            rejectVcf=$flareOutPath/$sample.chr$i.anc.${hapName}_unmapped.vcf.gz

            gatk LiftoverVcf \
                --CHAIN $chain \
                --INPUT $flareVcf \
                --OUTPUT $liftedVcf \
                --REFERENCE_SEQUENCE $referenceFasta \
                --REJECT $rejectVcf \
                --RECOVER_SWAPPED_REF_ALT true \
                &> $flareOutPath/liftover.$sample.chr$i.$hapName.log
        done
    done
done
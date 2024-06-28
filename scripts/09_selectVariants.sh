#!/bin/bash

## GATK extract variants - SNPs
# input : rawVariants.vcf.gz : extracted raw variants (with snps and indels together)
# output : vcf with only SNPs or INDELs

#SBATCH --job-name=selectSNPs
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=selectSNPs_INDELs.out
#SBATCH --error=selectSNPs_INDELs.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE

source ~/.bashrc
module load openjdk/11.0.2

export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna
export ANAL=$PROJECT/09_selectVariants

cd $ANAL
export INPUT=$PROJECT/08_genotypeGVCFs/rawVariants.vcf.gz
export OUTPUT=$ANAL/cod3_rawSNPs.vcf
## SNPs
java -Xmx96G -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
      -R $REFGEN \
      -V $INPUT \
      --select-type-to-include SNP \
      -O $OUTPUT

# INDELs
export OUTPUT=$ANAL/cod3_rawINDELs.vcf
java -Xmx96G -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
           -R $REFGEN \
           -V $INPUT \
           --select-type-to-include INDEL \
           -O $OUTPUT

jobinfo

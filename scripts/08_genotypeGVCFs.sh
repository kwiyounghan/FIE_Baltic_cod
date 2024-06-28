#!/bin/bash

## GATK GenotypeGVCFs  - joint genotype calling
# input : gvcf work directory from genomicsdbimport
# output : A final VCF in which all samples have been jointly genotyped.

#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=8-00:00:00
#SBATCH --output=GenotypeGVCFs.out
#SBATCH --error=GenotypeGVCFs.err
#SBATCH --partition=cluster
#SBATCH --qos=long
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=ALL
#SBATCH --export=NONE


source ~/.bashrc
module load openjdk/11.0.2

export WORK2=/gxfs_work2/geomar/smomw426
export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export GATK=$HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar
export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna
export ANAL=$PROJECT/08_genotypeGVCFs

cd $PROJECT/08_genotypeGVCFs/

java -Xmx96G -jar $GATK GenotypeGVCFs \
   -R $REFGEN --tmp-dir $WORK2/temp_dir \
   -V gendb://$PROJECT/07_genomicsdbimport/genomicdbimport\
   -O rawVariants.vcf.gz

jobinfo

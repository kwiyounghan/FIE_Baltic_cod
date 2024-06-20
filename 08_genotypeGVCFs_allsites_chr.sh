#!/bin/bash

## GATK GenotypeGVCFs  - joint genotype calling
# input : gvcf work directory from genomicsdbimport
# output : A final VCF in which all samples have been jointly genotyped.

#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=10-00:00:00
#SBATCH --output=GenotypeGVCFs_allsites.out
#SBATCH --error=GenotypeGVCFs_allsites.err
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
export ANAL=$PROJECT/08_genotypeGVCFs_chr

export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna


## since runing genotypeGVCFs with whole genome at once takes too long, break the jobs chr by chr
## chromosome information in $PROJECT/ref/intervals_allchr_no_unplaced.list) each line with a chromosome name from the reference genome
cd $ANAL

for line in $(cat $PROJECT/ref/intervals_allchr_no_unplaced.list); do
  java -Xmx400G -jar $GATK GenotypeGVCFs \
  -R $REFGEN --tmp-dir $WORK2/temp_dir -all-sites TRUE \
  -V gendb://$PROJECT/07_genomicsdbimport_chr/genomicdbimport_${line} \
  -O all_callable_sites_${line}.vcf.gz -L ${line} --genomicsdb-shared-posixfs-optimizations true

done

jobinfo

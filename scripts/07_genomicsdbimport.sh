#!/bin/bash

## GATK GenomicsDBImport : conslidating gvcf files

#SBATCH --job-name=genomicdbimport
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=10-00:00:00
#SBATCH --output=genomicdbimport.out
#SBATCH --error=genomicdbimport.err
#SBATCH --partition=cluster
#SBATCH --qos=long
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE

source ~/.bashrc
module load openjdk/11.0.2

export WORK2=/gxfs_work2/geomar/smomw426
export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export GATK=$HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar
export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna
export INBAM=$PROJECT/bams_233 #this is a list of _clean_dedup.bam files with absolute paths
export ANAL=$PROJECT/07_genomicsdbimport
#sample.map made by a script called make_sample_map.sh
#cd $ANAL
#cp $WORK/cod_all_5_115/07_genomicsdbimport/intervals_allchr.list .
cd $ANAL


    java -Xmx128G -jar $GATK GenomicsDBImport \
        --sample-name-map $ANAL/sample.map \
        --genomicsdb-workspace-path $ANAL/genomicdbimport/ \
        --tmp-dir $WORK2/temp_dir -L $PROJECT/ref/intervals_allchr.list --genomicsdb-shared-posixfs-optimizations true


jobinfo

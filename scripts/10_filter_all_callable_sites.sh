#!/bin/bash

#SBATCH --job-name=cod_invariant_filter
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=6-00:00:00
#SBATCH --output=invariant_filter_for_pi.out
#SBATCH --error=invariant_filter_for_pi.err
#SBATCH --partition=cluster
#SBATCH --qos=long
## qos - normal, long (10days), express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE



source ~/.bashrc
export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export ANAL=$PROJECT/09_selectVariants
export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna

cd $ANAL
#for GATK
module load openjdk/11.0.2

#https://pixy.readthedocs.io/en/latest/guide/pixy_guide.html#run-pixy
#for all sites vcf, still needs filtering of sites based on the quality etc.

export INPUT=$PROJECT/08_genotypeGVCFs/all_callable_sites.vcf.gz
export OUTPUT=$ANAL/all_callable_sites_flagged.vcf

java -Xmx48G -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration \
   -R $REFGEN \
   -V $INPUT \
   --filter-expression "QD < 4.0" \
   --filter-name "low_QD" \
   --filter-expression "QUAL <  30.0" \
   --filter-name "low_QUAL" \
   --filter-expression "FS > 25.0" \
   --filter-name "high_FS" \
   --filter-expression "SOR > 2.0" \
   --filter-name "high_SOR" \
   --filter-expression "MQ < 55.0" \
   --filter-name "low_MQ" \
   --filter-expression "MQRankSum < -1.5" \
   --filter-name "low_MQRankSum" \
   --filter-expression "ReadPosRankSum < -2.0" \
   --filter-name "low_ReadPosRankSum" \
   -O $OUTPUT


   # ------  kick out unfiltered variants -------------

       export INPUT=$ANAL/all_callable_sites_flagged.vcf
       export OUTPUT=$ANAL/all_callable_sites_filtered.vcf

       java -Xmx48G -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
           -R $REFGEN \
           -V $INPUT \
           -O $OUTPUT \
           --exclude-filtered


rm $INPUT
# create a filtered VCF containing only invariant sites
module load miniconda3
conda activate pop #for vcftools

export VCF=$ANAL/all_callable_sites_filtered.vcf
export VCFINV=all_callable_sites_filtered_invariant_miss0.8_dp10.vcf.gz
vcftools --gzvcf $VCF \
--max-missing 0.8 \
--min-meanDP 10 \
--max-meanDP 500 \
--max-maf 0 \
--recode --recode-INFO-all --stdout | bgzip -c > $VCFINV

tabix $VCFVAR
rm $VCF


jobinfo

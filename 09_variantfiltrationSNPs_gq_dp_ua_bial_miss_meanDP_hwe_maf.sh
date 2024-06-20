#!/bin/bash

##filtering snps with vcftools

#SBATCH --job-name=variantfiltration
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=variantfiltration_vcftools_SNPs_gq20_dp10_ua_mis08_bial_hwe_meanDP_maf.out
#SBATCH --error=variantfiltration_vcftools_SNPs_gq20_dp10_ua_mis08_bial_hwe_meanDP_maf.err
#SBATCH --partition=cluster
#SBATCH --qos=normal
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
   #SNPs

module load miniconda3
conda activate pop

export INPUT=$ANAL/cod3_filteredSNPs.vcf

#filter for quality --> genotype quality 20
export INPUT=$ANAL/cod3_filteredSNPs.vcf
export OUTPUT="cod3_filteredSNPs_gq20"
export FILTER="--minGQ 20"
vcftools --vcf $INPUT --recode --recode-INFO-all --out $OUTPUT $FILTER

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT

#filter for depth of 8
export INPUT=$ANAL/${OUTPUT}.recode.vcf
export OUTPUT="cod3_filteredSNPs_gq20_dp8"
export FILTER="--minDP 8"
vcftools --vcf $INPUT --recode --recode-INFO-all --out $OUTPUT $FILTER

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT

#filter unused alternates
export INPUT=$ANAL/${OUTPUT}.recode.vcf
export OUTPUT=${OUTPUT}_ua
module load openjdk/11.0.2
java -Xmx48g -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
--remove-unused-alternates TRUE \
--exclude-non-variants TRUE \
-R $REFGEN \
-V $INPUT \
-O ${ANAL}/${OUTPUT}.vcf

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT

#biallelic –min-alleles 2 –max-alleles 2
# this is possible because we excluded sites with "removeUnusedAlternates"
export INPUT=$ANAL/${OUTPUT}.vcf
export OUTPUT=${OUTPUT}_bial
vcftools --vcf $INPUT --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --out $OUTPUT

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT

#filter on missingness
export INPUT=$ANAL/${OUTPUT}.recode.vcf
export OUTPUT=${OUTPUT}_miss08
export FILTER="--max-missing 0.8"
vcftools --vcf $INPUT --recode --recode-INFO-all --remove-filtered-all --out $OUTPUT $FILTER

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT

#filter on minimum and maximum mean Depth for sites
export INPUT=$ANAL/${OUTPUT}.recode.vcf
export OUTPUT=${OUTPUT}_meanDP10
export FILTER="--min-meanDP 10"
vcftools --vcf $INPUT --recode --recode-INFO-all --remove-filtered-all --out $OUTPUT $FILTER

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT


#filter on excess of heterozygosity p value < 0.001
export INPUT=$ANAL/${OUTPUT}.recode.vcf
export OUTPUT=${OUTPUT}_hwe001
export FILTER="--hwe 0.001"
vcftools --vcf $INPUT --recode --recode-INFO-all --out $OUTPUT $FILTER

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT

#filter for MAF 5% = 0.05
export INPUT=$ANAL/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001.recode.vcf
export FILTER="--maf 0.05"
vcftools --vcf $INPUT --recode --recode-INFO-all --out ${OUTPUT}_maf0.05 $FILTER
#filter for MAF 5% = 0.01
export FILTER="--maf 0.01"
vcftools --vcf $INPUT --recode --recode-INFO-all --out ${OUTPUT}_maf0.01 $FILTER
#filter for MAF 5% = 0.005
export FILTER="--maf 0.005"
vcftools --vcf $INPUT --recode --recode-INFO-all --out ${OUTPUT}_maf0.005 $FILTER

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT

echo $INPUT
grep -v "^#" $INPUT | wc -l
rm $INPUT


jobinfo

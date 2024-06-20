#!/bin/bash

#Subsetting the master vcf file with 231 individuals
#SBATCH --job-name=cod_subset_vcf
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=cod_subset_vcf.out
#SBATCH --error=cod_subset_vcf.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long (10days), express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE



source ~/.bashrc
export WORK2=/gxfs_work2/geomar/smomw426
export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export ANAL=$PROJECT/09_selectVariants
export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna

cd $ANAL
module load miniconda3; conda activate sambcfenv

export VCF=$ANAL/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005.recode.vcf
# This vcf file includes 231 individuals from my study temporal and phenotypic individuals of EBC (115+47=162) + BOR (23),KIE (22),NOR(24) total 69indi. from Barth et al. 2019
# subset1 temporal EBC (115) +BOR(23)+KIE(22) (160)
# subset2 temporal EBC only (115) #excluding 1 individual from temporal EBC that I know is WBC
# subset3 temporal (115) + phenotypic EBC (47) =162

#groups_file.txt was created in my local computer, in metadata folder groups_file.txt
bcftools +split $VCF -Ov --groups-file $WORK2/cod3/09_selectVariants/groups_file.txt -o $ANAL

#then remove fixed sites
#"EBC_WBC_160.vcf","tempEBC_115.vcf","allEBC_164.vcf"
export INPUT=EBC_WBC_160.vcf
export OUTPUT=cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005_ebc_wbc_160/EBC_WBC_160_clean.vcf
module load openjdk/11.0.2
java -Xmx48g -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
--remove-unused-alternates TRUE \
--exclude-non-variants TRUE \
-R $REFGEN \
-V $INPUT \
-O $OUTPUT

export INPUT=tempEBC_115.vcf
export OUTPUT=cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005_temp115/tempEBC_115_clean.vcf
java -Xmx48g -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
--remove-unused-alternates TRUE \
--exclude-non-variants TRUE \
-R $REFGEN \
-V $INPUT \
-O $OUTPUT

export INPUT=allEBC_164.vcf
export OUTPUT=cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005_allEBC_164/allEBC_164_clean.vcf
java -Xmx48g -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
--remove-unused-alternates TRUE \
--exclude-non-variants TRUE \
-R $REFGEN \
-V $INPUT \
-O $OUTPUT


jobinfo

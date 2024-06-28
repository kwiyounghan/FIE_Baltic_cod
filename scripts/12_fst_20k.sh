#!/bin/bash

##project cod3
## performs Fst analysis - windowed in vcftools
##the main commands created by make_Fst.sh to make all population pairs

#SBATCH --job-name=Fst
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=cod3_Fst_window_20k.out
#SBATCH --error=cod3_Fst_window_20k.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE


source ~/.bashrc
export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export ANAL=$PROJECT/13_fst
cd $ANAL
#for vcftools to work.. install vcftools in a conda env
module load miniconda3
conda activate pop

export vcf=/gxfs_work2/geomar/smomw426/cod3/09_selectVariants/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005_temp115/tempEBC_115_clean.vcf

## windowed
#paired
export pop1=/gxfs_work2/geomar/smomw426/cod3/09_selectVariants/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005_temp115/1996.list
export pop2=/gxfs_work2/geomar/smomw426/cod3/09_selectVariants/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005_temp115/2019.list
export outpre="persite_fst_1996_2019_20k"
vcftools --gzvcf $vcf --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out $outpre  --fst-window-size 20000 --fst-window-step 20000


jobinfo

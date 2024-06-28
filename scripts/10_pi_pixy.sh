#!/bin/bash

#SBATCH --job-name=pi_pixy
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=pi_pixy.out
#SBATCH --error=pi_pixy.err
#SBATCH --partition=base
#SBATCH --qos=normal
## qos - normal, long (10days), express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE



source ~/.bashrc
export PROJECT=/gxfs_work/geomar/smomw426/cod3
export ANAL=$PROJECT/

module load gcc12-env/12.3.0
module load miniconda3/23.5.2

conda activate pop
cd $PROJECT/09_selectVariants/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005
export VCF=$PROJECT/09_selectVariants/cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005/invariant_cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.005.vcf.gz
tabix $VCF
conda deactivate

cd $ANAL
conda activate pixyenv

## newly installing pixy in a conda environment gave me an error "AttributeError: module 'numpy' has no attribute 'typeDict'"
## installing a specific numpy version in the env solved the error
# pip install numpy==1.22.1

## calculate pi
#this is needed for pixy program
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

#[pixy] pixy 1.2.7.beta1
export outtxt=invariant_cod3_filteredSNPs_gq20_dp8_ua_bial_miss08_meanDP10_hwe001_maf0.05
pixy --stats pi dxy --vcf $VCF --population $PROJECT/pop_info/pop_info --window_size 50000 --n_cores 32 --output_prefix $outtxt

jobinfo

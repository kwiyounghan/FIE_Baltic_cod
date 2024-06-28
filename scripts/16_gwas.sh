#!/bin/bash

#SBATCH --job-name=gwas
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=gwas_univar.out
#SBATCH --error=gwas_univar.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long (10days), express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE



source ~/.bashrc
export WORK2=/gxfs_work2/geomar/smomw426
export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export ANAL=$PROJECT/16_gwas

module load miniconda3
#conda create --name gwas gemma
conda activate gwas
cd $ANAL

## run with phi
# paste input/allEBC_152_clean.fam phenotype_linf_k_phi_sex.temp.csv | awk '{print $1,$2,$3,$4,$5,$11}' > allEBC_152_clean.fam
export bFile=/gxfs_work2/geomar/smomw426/cod3/16_gwas/input_uni_phi/allEBC_152_clean
export outpre=maf0.005_allEBC_152_clean_phi_temp
gemma -bfile $bFile -k output0/maf0.005_allEBC_152_clean_linfi_temp.cXX.txt -lmm 4 -o $outpre -c covariates.txt -miss 0.1 -maf 0.05


## from log file "maf0.005_allEBC_152_clean_linfi_temp.log.txt"
## ## GEMMA Version    = 0.98.3 (2020-11-28)
## number of total SNPs/var = 6273417
## number of analyzed SNPs/var = 679584
## REMLE log-likelihood in the null model = -270.396
## MLE log-likelihood in the null model = -268.944
## pve estimate in the null model = 0.195226
## se(pve) in the null model = 0.641648
## vg estimate in the null model = 3.462
## ve estimate in the null model = 1.7336
## beta estimate in the null model =   7.66888  0.302776
## se(beta) =   0.167128  0.238297

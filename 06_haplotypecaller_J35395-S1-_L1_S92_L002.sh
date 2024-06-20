#!/bin/bash
#SBATCH --job-name=haplotypecaller_J35395-S1-_L1_S92_L002
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=36:00:00
#SBATCH --output=haplotypecaller_J35395-S1-_L1_S92_L002.out
#SBATCH --error=haplotypecaller_J35395-S1-_L1_S92_L002.err
#SBATCH --partition=cluster
#SBATCH --qos=normal
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE


cd
source ~/.bashrc
module load openjdk/11.0.2


java -Xmx96G -jar /gxfs_home/geomar/smomw426/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar HaplotypeCaller \
  -R /gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna \
  -I /gxfs_work2/geomar/smomw426/cod3/04_dedup/J35395-S1-_L1_S92_L002_clean_dedup.bam \
  -O /gxfs_work2/geomar/smomw426/J35395-S1-_L1_S92_L002_output.raw.snps.indels.g.vcf -ERC GVCF --tmp-dir /gxfs_work2/geomar/smomw426/temp_dir

jobinfo

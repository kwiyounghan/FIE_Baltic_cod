#!/bin/bash
#SBATCH --job-name=MarkIlluminaAdapters
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=markilluminaadapters_J35395-S1-L1_S92.out
#SBATCH --error=markilluminaadapters_J35395-S1-L1_S92.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

cd /gxfs_work1/geomar/smomw426/cod_mod2/02_markilluminaadapters
source ~/.bashrc
module load openjdk/11.0.2

echo /gxfs_work1/geomar/smomw426/cod_mod2/01_fastq2ubam/J35395-S1-L1_S92_unmapped.bam marking adapters
java -Xmx72G -jar /gxfs_home/geomar/smomw426/software/picard.jar MarkIlluminaAdapters \
	I=/gxfs_work1/geomar/smomw426/cod_mod2/01_fastq2ubam/J35395-S1-L1_S92_unmapped.bam \
	O=/gxfs_work1/geomar/smomw426/cod_mod2/02_markilluminaadapters/J35395-S1-L1_S92_unmapped_xt.bam \
	M=/gxfs_work1/geomar/smomw426/cod_mod2/02_markilluminaadapters/J35395-S1-L1_S92_unmapped_xt_metrics.txt \
	TMP_DIR=/gxfs_work1/geomar/smomw426/temp_dir/

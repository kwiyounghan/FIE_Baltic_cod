#!/bin/bash
#SBATCH --job-name=dedup_J35395-S1-_L1_S92_L002
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=dedup_J35395-S1-_L1_S92_L002.out
#SBATCH --error=dedup_J35395-S1-_L1_S92_L002.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE

cd /gxfs_work2/geomar/smomw426/cod3/04_dedup
source ~/.bashrc
module load openjdk/11.0.2

echo /gxfs_work2/geomar/smomw426/cod3/03_map_clean/J35395-S1-L1_S92_L002_clean.bam deduplicating
java -Xmx72G -jar /gxfs_home/geomar/smomw426/software/picard.jar MarkDuplicates \
      I=/gxfs_work2/geomar/smomw426/cod3/03_map_clean/J35395-S1-L1_S92_L002_clean.bam \
      O=J35395-S1-_L1_S92_L002_clean_dedup.bam \
      CREATE_INDEX=true \
      METRICS_FILE=J35395-S1-_L1_S92_L002_clean_dedup_metrics.txt \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=/gxfs_work2/geomar/smomw426/temp_dir
jobinfo

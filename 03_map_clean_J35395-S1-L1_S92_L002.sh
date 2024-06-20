#!/bin/bash
#SBATCH --job-name=map_clean_ones_J35395-S1-L1_S92_L002
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=map_clean_J35395-S1-L1_S92_L002.out
#SBATCH --error=map_clean_J35395-S1-L1_S92_L002.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE
cd /gxfs_work1/geomar/smomw426/cod3/03_map_clean
set -o pipefail
source ~/.bashrc
module load openjdk/11.0.2

java -Xmx72g -jar /gxfs_home/geomar/smomw426/software/picard.jar SamToFastq \
INPUT=/gxfs_work1/geomar/smomw426/cod_three_and_barth/02_markilluminaadapters/J35395-S1-L1_S92_L002_unmapped_xt.bam FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true \
TMP_DIR=/gxfs_work2/geomar/smomw426/temp_dir | \

/gxfs_home/geomar/smomw426/software/bwa/bwa mem -M -t 32 -p /gxfs_work1/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna /dev/stdin | \

java -Xmx72G -jar /gxfs_home/geomar/smomw426/software/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=/gxfs_work1/geomar/smomw426/cod_three_and_barth/01_fastq2ubam/J35395-S1-L1_S92_L002_unmapped.bam \
OUTPUT= /gxfs_work1/geomar/smomw426/cod3/03_map_clean/J35395-S1-L1_S92_L002_clean.bam \
R=/gxfs_work1/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna \
MAX_INSERTIONS_OR_DELETIONS=-1 CREATE_INDEX=true ATTRIBUTES_TO_RETAIN=XS PRIMARY_ALIGNMENT_STRATEGY=MostDistant CLIP_ADAPTERS=false ADD_MATE_CIGAR=true TMP_DIR=/gxfs_work2/geomar/smomw426/temp_dir
echo done

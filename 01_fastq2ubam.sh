#!/bin/bash

## 01_fastq2ubam

#SBATCH --job-name=fastq2ubam
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=fastq2ubam.out
#SBATCH --error=fastq2ubam.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

cd /gxfs_work1/geomar/smomw426/cod_mod2/01_fastq2ubam
source ~/.bashrc
module load openjdk/11.0.2
export PROJECT=/gxfs_work1/geomar/smomw426/cod_mod2
export RAW_FQ=/gxfs_work1/geomar/smomw426/cod_mod2/raw/*.fastq.gz
#export REFGEN=

for i in $RAW_FQ; do
	# (F12149-L1)_(S1)_(L001)_R1_001.fastq.gz
	#F12149-L1 bashrematch1 S1 bashrematch2 L001 bashrematch3
  # J35433-S1-L1_S39_L001_R1_001.fastq
  # J35433-S1-L1_S39_L001_R2_001.fastq
     if [[ $i =~ /gxfs_work1/geomar/smomw426/cod_mod2/raw/(.*)_(.*)_(.*)_R1_001.fastq.gz ]];
     then
     echo $i test1 ${BASH_REMATCH[1]}_${BASH_REMATCH[2]} test2 ${BASH_REMATCH[1]} bashrematch1 ${BASH_REMATCH[2]} bashrematch2 ${BASH_REMATCH[3]} bashrematch3

     java -Xmx72G -jar /gxfs_home/geomar/smomw426/software/picard.jar FastqToSam \
	FASTQ= $i \
	FASTQ2= $PROJECT/raw/${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_${BASH_REMATCH[3]}_R2_001.fastq.gz \
	OUTPUT= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_unmapped.bam \
	READ_GROUP_NAME= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]} \
	SAMPLE_NAME= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]} \
	LIBRARY_NAME= ${BASH_REMATCH[1]}_${BASH_REMATCH[2]} \
	PLATFORM_UNIT= ${BASH_REMATCH[3]} \
	PLATFORM= ILLUMINA TMP_DIR=$WORK/temp_dir/

  fi

done

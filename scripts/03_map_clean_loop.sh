#!/bin/bash
#this script creates files for each individual, one unmapped_xt bam file, and runs those files after creating them.
#piping through picard SamToFastq - bwa - picard MergeBamAlignment to creat "clean" bam files
#input is unmapped_xt.bam and output is _clean.bam
#mkdir $PROJECT/scripts/03_map_clean_ones
#readmapping for my cod seq ALL individuals 1996-1998, 2002, 2008, 2014, 2019 with Gadmor3 ref

## index the reference file before readmapping!  /gxfs_home/geomar/smomw426/software/bwa/bwa index $REFGEN
##seq dictionary present?  /gxfs_home/geomar/smomw426/software/picard.jar CreateSequenceDictionary -R $REFGEN
#run in the working directory bash 03_map_clean_one_loop.sh , not sbatch.. .
#make a script that checks all the prerequisite for the runs?!

export PROJECT=/gxfs_work1/geomar/smomw426/cod3
export ANAL=$PROJECT/03_map_clean
cd $ANAL

#reference genome
export REFGEN=/gxfs_work1/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna
#/gxfs_home/geomar/smomw426/software/bwa/bwa index $REFGEN
#java -jar /gxfs_home/geomar/smomw426/software/picard.jar CreateSequenceDictionary -R $REFGEN


########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
#samples from NSC 1996, 1996-1998, 2002, 2008
export PROJECT_ub=/gxfs_work1/geomar/smomw426/cod_1996_2002_2008_and_phenotype_samples
export uBAM_XT=$PROJECT_ub/02_markilluminaadapters/*_xt.bam
#108 uBAM_XT 1996, 1996-1998, 2002, 2008
for i in $uBAM_XT ; do
    if [[ $i =~ $PROJECT_ub/02_markilluminaadapters/(.*)_unmapped_xt.bam ]]; # ${BASH_REMATCH[1]} J32532-L1_S1
    then
      if [[ ! -f "$ANAL/${BASH_REMATCH[1]}_clean.bam" ]]; then

		#create a file with shebang and add slurm info
		echo '#!/bin/bash' > $PROJECT/scripts/03_map_clean_ones/03_map_clean_${BASH_REMATCH[1]}.sh
		export file=$PROJECT/scripts/03_map_clean_ones/03_map_clean_${BASH_REMATCH[1]}.sh
    echo "#SBATCH --job-name=map_clean_ones_${BASH_REMATCH[1]}"  >> $file
    echo '#SBATCH --nodes=1'  >> $file
    echo '#SBATCH --tasks-per-node=1'  >> $file
    echo '#SBATCH --cpus-per-task=32'  >> $file
    echo '#SBATCH --mem=96G'  >> $file
    echo '#SBATCH --time=12:00:00'  >> $file
    echo "#SBATCH --output=map_clean_${BASH_REMATCH[1]}.out"  >> $file
    echo "#SBATCH --error=map_clean_${BASH_REMATCH[1]}.err"  >> $file
    echo '#SBATCH --partition=cluster'  >> $file
    echo '#SBATCH --qos=express'  >> $file
    echo '## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos'  >> $file
    echo '#SBATCH --mail-user=khan@geomar.de'  >> $file
    echo '#SBATCH --mail-type=FAIL'  >> $file
    echo '#SBATCH --export=NONE'  >> $file

    #load any software you need !
		#add commands and options
		#set working directory and environment
		echo "cd $ANAL" >> $file
		echo 'set -o pipefail' >> $file
		echo 'source ~/.bashrc' >> $file
    echo 'module load openjdk/11.0.2' >> $file
    echo '' >> $file

		#SamToFastq
		echo "java -Xmx72g -jar /gxfs_home/geomar/smomw426/software/picard.jar SamToFastq \\" >> $file
		echo "INPUT=$i FASTQ=/dev/stdout \\" >> $file
		echo "CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true \\" >> $file
		echo "TMP_DIR=$WORK2/temp_dir | \\" >> $file
    echo '' >> $file
		#bwa in $HOME/software/
		echo "/gxfs_home/geomar/smomw426/software/bwa/bwa mem -M -t 32 -p $REFGEN /dev/stdin | \\" >> $file
    echo '' >> $file
		#MergeBamAlignment
		echo "java -Xmx72G -jar /gxfs_home/geomar/smomw426/software/picard.jar MergeBamAlignment \\" >> $file
		echo "ALIGNED_BAM=/dev/stdin \\" >> $file
		echo "UNMAPPED_BAM=$PROJECT_ub/01_fastq2ubam/${BASH_REMATCH[1]}_unmapped.bam \\" >> $file
		echo "OUTPUT= $PROJECT/03_map_clean/${BASH_REMATCH[1]}_clean.bam \\" >> $file
		echo "R=$REFGEN \\" >> $file
		echo "MAX_INSERTIONS_OR_DELETIONS=-1 CREATE_INDEX=true ATTRIBUTES_TO_RETAIN=XS PRIMARY_ALIGNMENT_STRATEGY=MostDistant CLIP_ADAPTERS=false ADD_MATE_CIGAR=true TMP_DIR=$WORK2/temp_dir" >> $file
		echo "echo done" >> $file

		sbatch $file
    sleep 3s


    fi
	fi
done

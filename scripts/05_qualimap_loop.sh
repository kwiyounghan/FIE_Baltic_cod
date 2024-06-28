#!/bin/bash
##this script creates a file for each individual, , and runs those files after creating them.
# for qualimap : input is inbam files
# mkdir $PROJECT/scripts/05_qualimap_ones
#check if all the software is installed.. in conda activate pop

export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export INBAM=$PROJECT/04_dedup/*_clean_dedup.bam

for i in $INBAM ; do
    if [[ $i =~ $PROJECT/04_dedup/(.*)_clean_dedup.bam ]]; then # ${BASH_REMATCH[1]} J32532-L1_S1
    if [[ ! -d "$PROJECT/05_QC/qualimap/${BASH_REMATCH[1]}" ]]; then

      echo '#!/bin/bash' > $PROJECT/scripts/05_qualimap_ones/${BASH_REMATCH[1]}.sh
      export file=$PROJECT/scripts/05_qualimap_ones/${BASH_REMATCH[1]}.sh
      echo "" >> $file
      echo "#This script uses qualimap in conda." >> $file

      echo "#SBATCH --job-name=05_qualimap_${BASH_REMATCH[1]}"  >> $file
      echo '#SBATCH --nodes=1'  >> $file
      echo '#SBATCH --tasks-per-node=1'  >> $file
      echo '#SBATCH --cpus-per-task=32'  >> $file
      echo '#SBATCH --mem=16G'  >> $file
      echo '#SBATCH --time=4:00:00'  >> $file
      echo "#SBATCH --output=$PROJECT/05_QC/qualimap/${BASH_REMATCH[1]}.out"  >> $file
      echo "#SBATCH --error=$PROJECT/05_QC/qualimap/${BASH_REMATCH[1]}.err"  >> $file
      echo '#SBATCH --partition=cluster'  >> $file
      echo '#SBATCH --qos=express'  >> $file
      echo '## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos'  >> $file
      echo '#SBATCH --mail-user=khan@geomar.de'  >> $file
      echo '#SBATCH --mail-type=FAIL'  >> $file
      echo '#SBATCH --export=NONE'  >> $file
      echo '' >> $file


      echo "source ~/.bashrc" >> $file
      echo "cd $PROJECT/05_QC/qualimap" >> $file
      echo "module load miniconda3" >> $file
      echo "conda activate pop" >> $file
      echo "qualimap bamqc -bam $i -outdir $PROJECT/05_QC/qualimap/${BASH_REMATCH[1]} --java-mem-size=16G" >> $file


      sleep 2s
      sbatch $file
    fi
  fi
  done 

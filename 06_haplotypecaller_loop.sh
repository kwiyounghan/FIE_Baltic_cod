#!/bin/bash

#this script creates a file for each individual, , and runs those files after creating them.
# for HaplotypeCaller GATK
# the input file is _clean_dedup.bam and outbut files are _output.raw.snps.indels.g.vcf
#is reference file indexed ? .fai file exists?
# samtools faidx ref.fasta

export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export ANAL=$PROEJCT/06_haplotypecaller
export GATK=$HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar
export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna
export INBAM=$PROJECT/clean_dedup_bams_233 #this is a list of _clean_dedup.bam files with absolute paths

mkdir $PROJECT/scripts/06_haplotypecaller_ones
for i in $(cat $INBAM) ; do
    if [[ $i =~ (.*)/04_dedup/(.*)_clean_dedup.bam ]]; # ${BASH_REMATCH[1]} J32532-L1_S1
    then
      #slurm options
      echo '#!/bin/bash' > $PROJECT/scripts/06_haplotypecaller_ones/${BASH_REMATCH[2]}.sh
  		export file=$PROJECT/scripts/06_haplotypecaller_ones/${BASH_REMATCH[2]}.sh

      echo "#SBATCH --job-name=haplotypecaller_${BASH_REMATCH[2]}"  >> $file
      echo '#SBATCH --nodes=1'  >> $file
      echo '#SBATCH --tasks-per-node=1'  >> $file
      echo '#SBATCH --cpus-per-task=32'  >> $file
      echo '#SBATCH --mem=96G'  >> $file
      echo '#SBATCH --time=36:00:00'  >> $file
      echo "#SBATCH --output=haplotypecaller_${BASH_REMATCH[2]}.out"  >> $file
      echo "#SBATCH --error=haplotypecaller_${BASH_REMATCH[2]}.err"  >> $file
      echo '#SBATCH --partition=cluster'  >> $file
      echo '#SBATCH --qos=normal'  >> $file
      echo '## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos'  >> $file
      echo '#SBATCH --mail-user=khan@geomar.de'  >> $file
      echo '#SBATCH --mail-type=FAIL'  >> $file
      echo '#SBATCH --export=NONE'  >> $file
      echo '' >> $file
      echo '' >> $file

      # GATK HaplotypeCaller
      #load any software you need !
  		#add commands and options
  		#set working directory and environment
  		echo "cd $ANAL" >> $file
  		echo 'source ~/.bashrc' >> $file
      echo 'module load openjdk/11.0.2' >> $file
      echo '' >> $file
      echo '' >> $file
      echo "java -Xmx96G -jar $GATK HaplotypeCaller \\" >> $file
      echo "  -R $REFGEN \\" >> $file
      echo "  -I $i \\" >> $file
      echo "  -O /gxfs_work2/geomar/smomw426/${BASH_REMATCH[2]}_output.raw.snps.indels.g.vcf -ERC GVCF --tmp-dir /gxfs_work2/geomar/smomw426/temp_dir" >> $file
      echo '' >> $file
      echo 'jobinfo' >> $file

      sbatch $file
      sleep 2s
    fi
  done

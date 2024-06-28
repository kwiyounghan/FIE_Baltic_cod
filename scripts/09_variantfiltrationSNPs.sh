#!/bin/bash

##This creates a VCF file called filtered_snps.vcf, containing all the original SNPs from the rawSNPs.vcf file,
# but now the SNPs are annotated with either PASS or FILTER depending on whether or not they passed the filters.

# It has "same" amount of lines as unfiltered, probably with another header
# Next step to create VCF-tools compatible vcf files -> this cuts out FITER SNPs --> less lines!

# input : rawSNPs.vcf : rawSNP vcf file
# mid output  : a vcf with PASS/FILTER flags but contains all rawSNPs
# output : a vcf file with filtered SNPs only

#SBATCH --job-name=variantfiltration
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=variantfiltrationSNPs.out
#SBATCH --error=variantfiltrationSNPs.err
#SBATCH --partition=cluster
#SBATCH --qos=express
## qos - normal, long, express(max time 12h), test(max time 2h), to see the list >sacctmgr show qos
#SBATCH --mail-user=khan@geomar.de
#SBATCH --mail-type=FAIL
#SBATCH --export=NONE

source ~/.bashrc
module load openjdk/11.0.2

export PROJECT=/gxfs_work2/geomar/smomw426/cod3
export REFGEN=/gxfs_work2/geomar/smomw426/cod3/ref/GCF_902167405.1_gadMor3.0_genomic.fna
export ANAL=$PROJECT/09_selectVariants

cd $ANAL
   #SNPs

   export INPUT=$ANAL/cod3_rawSNPs.vcf
   export OUTPUT=$ANAL/cod3_flagged_rawSNPs.vcf

java -Xmx48G -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration \
   -R $REFGEN \
   -V $INPUT \
   --filter-expression "QD < 4.0" \
   --filter-name "low_QD" \
   --filter-expression "QUAL <  30.0" \
   --filter-name "low_QUAL" \
   --filter-expression "FS > 25.0" \
   --filter-name "high_FS" \
   --filter-expression "SOR > 2.0" \
   --filter-name "high_SOR" \
   --filter-expression "MQ < 55.0" \
   --filter-name "low_MQ" \
   --filter-expression "MQRankSum < -1.5" \
   --filter-name "low_MQRankSum" \
   --filter-expression "ReadPosRankSum < -2.0" \
   --filter-name "low_ReadPosRankSum" \
   -O $OUTPUT

   export INPUT=$ANAL/cod3_flagged_rawSNPs.vcf
   export OUTPUT=$ANAL/cod3_flagged_rawSNPs.table

   java -Xmx48G -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable \
   -R $REFGEN \
   -V $INPUT \
   -O $OUTPUT \
   -F ID -F CHROM -F POS -F FILTER -F QD -F FS \
   -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -F NO-CALL -F TRANSITION -F HET -F NCALLED -F InbreedingCoeff \
   --show-filtered


# ------  kick out unfiltered variants -------------

    export INPUT=$ANAL/cod3_flagged_rawSNPs.vcf
    export OUTPUT=$ANAL/cod3_filteredSNPs.vcf

    java -Xmx48G -jar $HOME/software/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
        -R $REFGEN \
        -V $INPUT \
        -O $OUTPUT \
        --exclude-filtered

rm $INPUT

jobinfo

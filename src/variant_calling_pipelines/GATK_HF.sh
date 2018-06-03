#!/bin/bash

###############################################################################
#      Script for performing SNP calling and hard filtering with GATK         #
#                                                                             #
#   More info at http://gatkforums.broadinstitute.org/gatk/discussion/2806/   #
###############################################################################

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then
    printf "Usage: $0 num_threads input_bam chromosome\n"
    exit -1
fi

num_threads=$1
printf "Number of threads: $num_threads\n"
input_bam=$2
printf "Input BAM file: $input_bam\n"
input_bai="$(printf $input_bam | sed 's/\.[^.]*$//')".bai
printf "Corresponding BAM index file: $input_bai\n"
log_txt="$input.bam.GATK_HF.log"
printf "Log file: $log_txt\n"
chromosome=$3
printf "Chromosome: $chromosome\n"

if [ -f $log_txt ]; then printf "Error: File $log_txt file already exists.\n"; exit -1; fi
if [ ! -f $input_bam ]; then printf "Error: Input BAM file $input_bam is not a regular file.\n"; exit -1; fi
if [ ! -f $input_bai ]; then printf "Error: BAM index file $input_bai is not a regular file.\n"; exit -1; fi

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/data/voges/muenteferi/GATK_bundle-2.8-b37"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_vcf="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_vcf="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
KG_vcf="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_vcf="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_vcf="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

if [ ! -f $ref_fasta ]; then printf "Error: File $ref_fasta is not a regular file.\n"; exit -1; fi
if [ ! -f $hapmap_vcf ]; then printf "Error: File $hapmap_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $omni_vcf ]; then printf "Error: File $omni_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $KG_vcf ]; then printf "Error: File $KG_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $dbsnps_vcf ]; then printf "Error: File $dbsnps_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $mills_vcf ]; then printf "Error: File $mills_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $indels_vcf ]; then printf "Error: File $indels_vcf is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
bowtie2="/project/dna/install/bowtie2-2.2.5/bowtie2"
java="/usr/bin/java"
java_opts=""
samtools="/project/dna/install/samtools-1.3/bin/samtools"

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
picard_jar="/project/dna/install/picard-tools-2.4.1/picard.jar"

if [ ! -x $bowtie2 ]; then printf "Error: Binary file $bowtie2 is not executable.\n"; exit -1; fi
if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
if [ ! -f $GenomeAnalysisTK_jar ]; then printf "Error: JAR file $GenomeAnalysisTK_jar is not a regular file.\n"; exit -1; fi
if [ ! -f $picard_jar ]; then printf "Error: JAR file $picard_jar is not a regular file.\n"; exit -1; fi

###############################################################################
#                         Variant calling with GATK                           #
###############################################################################

printf "[1/3] Variant calling\n"
$java $java_opts -jar $GenomeAnalysisTK_jar -T HaplotypeCaller -nct $num_threads -R $ref_fasta -L $chromosome -I $input_bam --dbsnp $dbsnps_vcf --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $input_bam.GATK.raw_variants.vcf &>>$log_txt

printf "[2/3] SNP extraction\n"
$java $java_opts -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_fasta -L $chromosome -V $input_bam.GATK.raw_variants.vcf -selectType SNP -o $input_bam.GATK.snps.vcf &>>$log_txt

printf "[3/3] Hard filtering\n"
filterExpression="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
filterName="GATK_Recommended"
$java $java_opts -jar $GenomeAnalysisTK_jar -T VariantFiltration -R $ref_fasta -L $chromosome -V $input_bam.GATK.snps.vcf --filterExpression "$filterExpression" --filterName "$filterName" -o $input_bam.GATK.snps.hard_filtered.vcf &>>$log_txt

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup ... "
#
printf "Done\n";


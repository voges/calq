#!/bin/bash

###############################################################################
#      Script for performing SNP calling and hard filtering with GATK         #
#                                                                             #
#   More info at http://gatkforums.broadinstitute.org/gatk/discussion/2806/   #
###############################################################################

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 4 ]; then
    printf "Usage: $0 num_threads log_txt input_bam chromosome\n"
    exit -1
fi

num_threads=$1
printf "Number of threads: $num_threads\n"
log_txt=$2
printf "Log file: $log_txt\n"
input_bam=$3
printf "Input BAM file: $input_bam\n"
chromosome=$4
printf "Chromosome: $chromosome\n"

printf "Checking log file $log_txt ... "
if [ -f $log_txt ]; then printf "log file already exists: $log_txt\n"; exit -1; fi
printf "OK\n"

printf "Checking input BAM file $input_bam ... "
if [ ! -f $input_bam ]; then printf "did not find input BAM file: $input_bam\n"; exit -1; fi
printf "OK\n"

input_bai="$(printf $input_bam | sed 's/\.[^.]*$//')".bai
printf "Checking corresponding BAM index file $input_bai ... "
if [ ! -f $input_bai ]; then printf "did not find BAM index file: $input_bai\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/data/voges/MPEG/GATK_bundle-2.8-b37"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_vcf="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_vcf="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
KG_vcf="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_vcf="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_vcf="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

printf "Checking GATK bundle ... "
if [ ! -f $ref_fasta ]; then printf "did not find $ref_fasta\n"; exit -1; fi
if [ ! -f $hapmap_vcf ]; then printf "did not find $hapmap_vcf\n"; exit -1; fi
if [ ! -f $omni_vcf ]; then printf "did not find $omni_vcf\n"; exit -1; fi
if [ ! -f $KG_vcf ]; then printf "did not find $KG_vcf\n"; exit -1; fi
if [ ! -f $dbsnps_vcf ]; then printf "did not find $dbsnps_vcf\n"; exit -1; fi
if [ ! -f $mills_vcf ]; then printf "did not find $mills_vcf\n"; exit -1; fi
if [ ! -f $indels_vcf ]; then printf "did not find $indels_vcf\n"; exit -1; fi
printf "OK\n"

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

printf "Checking executables ... "
if [ ! -x $bowtie2 ]; then printf "did not find $bowtie2\n"; exit -1; fi
if [ ! -x $java ]; then printf "did not find $java\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "did not find $samtools\n"; exit -1; fi
if [ ! -e $GenomeAnalysisTK_jar ]; then printf "did not find $GenomeAnalysisTK_jar\n"; exit -1; fi
if [ ! -e $picard_jar ]; then printf "did not find $picard_jar\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                         Variant calling with GATK                           #
###############################################################################

printf "[1/3] Variant calling ... "
$java $java_opts -jar $GenomeAnalysisTK_jar -T HaplotypeCaller -nct $num_threads -R $ref_fasta -L $chromosome -I $input_bam --dbsnp $dbsnps_vcf --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $input_bam.raw_variants.vcf &>>$log_txt
printf "OK\n"

printf "[2/3] SNP extraction ... "
$java $java_opts -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_fasta -L $chromosome -V $input_bam.raw_variants.vcf -selectType SNP -o $input_bam.snps.vcf &>>$log_txt
printf "OK\n"

printf "[3/3] Hard filtering ... "
filterExpression="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
filterName="GATK_Recommended"
$java $java_opts -jar $GenomeAnalysisTK_jar -T VariantFiltration -R $ref_fasta -L $chromosome -V $input_bam.snps.vcf --filterExpression "$filterExpression" --filterName "$filterName" -o $input_bam.snps.hard_filtered.vcf &>>$log_txt
printf "OK\n";

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


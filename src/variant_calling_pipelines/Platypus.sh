#!/bin/bash

###############################################################################
#            Script for performing variant calling with Platypus              #
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

if [ ! -f $ref_fasta ]; then printf "Error: File $ref_fasta is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
java="/usr/bin/java"
java_opts=""
python="/usr/bin/python"

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
Platypus_py="/project/dna/install/Platypus-0.8.1/Platypus.py"

if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -f $GenomeAnalysisTK_jar ]; then printf "Error: JAR file $GenomeAnalysisTK_jar is not a regular file.\n"; exit -1; fi
if [ ! -f $Platypus_py ]; then printf "Error: JAR file $Platypus_py is not a regular file.\n"; exit -1; fi

printf "Sourcing project_dna.config\n"
source /project/dna/project_dna.config

###############################################################################
#                         Variant calling with GATK                           #
###############################################################################

printf "[1/2] Variant calling\n"
$python $Platypus_py callVariants --nCPU=$num_threads --bamFiles=$input_bam --refFile=$ref_fasta --regions=$chromosome --output=$input_bam.raw_variants.vcf --logFileName=$log_txt &>>$log_txt

printf "[2/2] SNP extraction\n"
$java $java_opts -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_fasta -L $chromosome -V $input_bam.raw_variants.vcf -selectType SNP -o $input_bam.snps.vcf &>>$log_txt

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n";


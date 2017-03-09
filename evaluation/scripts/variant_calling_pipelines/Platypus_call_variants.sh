#!/bin/bash

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

gatk_bundle_path="/phys/intern2/tmp/data_gidb/MPEG/GATK_bundle-2.8-b37"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"

printf "Checking GATK bundle ... "
if [ ! -f $ref_fasta ]; then printf "did not find $ref_fasta\n"; exit -1; fi
if [ ! -f $hapmap_vcf ]; then printf "did not find $hapmap_vcf\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

python="/usr/bin/python"
Platypus_py="/project/dna/install/Platypus-0.8.1/Platypus.py"

printf "Checking executables ... "
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -e $Platypus_py ]; then printf "did not find $Platypus_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                         Variant calling with GATK                           #
###############################################################################

printf "[1/1] Variant calling ... "
$python $Platypus --nCPU=$num_threads --bamFiles=$input_bam --refFile=$ref_fasta --regions=$chromosome --output=$input_bam.Platypus_snps.vcf

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


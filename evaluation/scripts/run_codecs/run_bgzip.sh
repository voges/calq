#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then
    printf "Usage: $0 input_fastq num_threads\n"
    exit -1
fi

input_fastq=$1
printf "Input FASTQ file: $input_fastq\n"
num_threads=$2
printf "Number of threads: $num_threads\n"

printf "Checking input FASTQ file $input_sam ... "
if [ ! -f $input_fastq ]; then printf "did not find input FASTQ file: $input_fastq\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
bgzip="/project/dna/install/htslib-1.3/bin/bgzip"
python="/usr/bin/python"

# Python scripts
xtract_qual_fastq_py="/home/voges/git/ngstools/xtract_qual_fastq.py"

printf "Checking executables ... "
if [ ! -x $bgzip ]; then printf "did not find $bgzip\n"; exit -1; fi
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -e $xtract_qual_fastq_py ]; then printf "did not find $xtract_qual_fastq_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Extracting quality values from FASTQ file ... "
$python $xtract_qual_fastq_py $input_fastq 2> $input_fastq.qual
printf "OK\n"

printf "Running bgzip ... "
$bgzip -@ $num_threads -c $input_fastq.qual > $input_fastq.qual.bgz
wc -c $input_fastq.qual.bgz > $input_fastq.qual.bgz.log
printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


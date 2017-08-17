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

printf "Checking input FASTQ file $input_fastq ... "
if [ ! -f $input_fastq ]; then printf "did not find input FASTQ file: $input_fastq\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
dsrc="/project/dna/install/dsrc-2.0/dsrc"
dsrc_string="dsrc_il8b"
gzip="/usr/bin/gzip"
python="/usr/bin/python"

# Python scripts
xtract_part_fastq_py="/home/voges/git/calq/src/ngstools/xtract_part_fastq.py"

printf "Checking executables ... "
if [ ! -x $dsrc ]; then printf "did not find $dsrc\n"; exit -1; fi
if [ ! -x $gzip ]; then printf "did not find $gzip\n"; exit -1; fi
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -e $xtract_part_fastq_py ]; then printf "did not find $xtract_part_fastq_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Compressing with DSRC ... "
$dsrc c -d3 -q2 -b256 -l -t$num_threads $input_fastq $input_fastq.$dsrc_string
printf "OK\n"

printf "Decompressing with DSRC ... "
$dsrc d -t$num_threads $input_fastq.$dsrc_string $input_fastq.$dsrc_string.fastq
printf "OK\n"

printf "Extracting quality values from FASTQ file ... "
$python $xtract_part_fastq_py $input_fastq.$dsrc_string.fastq 3 1> $input_fastq.$dsrc_string.fastq.qual
printf "OK\n"

printf "Running gzip ... "
$gzip -9 -c $input_fastq.$dsrc_string.fastq.qual > $input_fastq.$dsrc_string.fastq.qual.gz
wc -c $input_fastq.$dsrc_string.fastq.qual.gz > $input_fastq.$dsrc_string.fastq.qual.gz.log
printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


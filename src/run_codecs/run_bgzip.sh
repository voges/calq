#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 input_fastq num_threads\n"; exit -1; fi

input_fastq=$1
printf "Input FASTQ file: $input_fastq\n"
num_threads=$2
printf "Number of threads: $num_threads\n"

if [ ! -f $input_fastq ]; then printf "Error: Input FASTQ file $input_fastq is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
bgzip="/project/dna/install/htslib-1.3/bin/bgzip"
python="/usr/bin/python"

# Python scripts
xtract_part_fastq_py="/home/voges/git/calq/src/ngstools/xtract_part_fastq.py"

if [ ! -x $bgzip ]; then printf "Error: Binary file $bgzip is not executable.\n"; exit -1; fi
if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -f $xtract_part_fastq_py ]; then printf "Error: Python script $xtract_part_fastq_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Extracting quality values\n  from: $input_fastq\n  to: $input_fastq.qual\n"
$python $xtract_part_fastq_py $input_fastq 3 1> $input_fastq.qual

printf "Compressing with bgzip\n  from: $input_fastq.qual\n  to:  $input_fastq.qual.bgz\n"
$bgzip -@ $num_threads -c $input_fastq.qual 1> $input_fastq.qual.bgz

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


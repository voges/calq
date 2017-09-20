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
dsrc="/project/dna/install/dsrc-2.0/dsrc"
dsrc_string="dsrc_il8b"
gzip="/usr/bin/gzip"
python="/usr/bin/python"
time="/usr/bin/time"

# Python scripts
xtract_part_fastq_py="/home/voges/git/calq/src/ngstools/xtract_part_fastq.py"

if [ ! -x $dsrc ]; then printf "Error: Binary file $dsrc is not executable.\n"; exit -1; fi
if [ ! -x $gzip ]; then printf "Error: Binary file $gzip is not executable.\n"; exit -1; fi
if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi
if [ ! -f $xtract_part_fastq_py ]; then printf "Error: Python script $xtract_part_fastq_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Compressing with DSRC\n  from: $input_fastq\n  to: $input_fastq.$dscr_string\n"
if [ -f $input_fastq.$dscr_string ]; then
    printf "$input_fastq.$dscr_string already exists (not reproducing it)\n"
else
    $dsrc c -d3 -q2 -b256 -l -t$num_threads $input_fastq $input_fastq.$dsrc_string
fi

printf "Decompressing with DSRC\n  from: $input_fastq.$dscr_string\n  to: $input_fastq.$dsrc_string.fastq\n"
if [ -f $input_fastq.$dscr_string.fastq ]; then
    printf "$input_fastq.$dscr_string.fastq already exists (not reproducing it)\n"
else
    $dsrc d -t$num_threads $input_fastq.$dsrc_string $input_fastq.$dsrc_string.fastq
fi

printf "Extracting quality values\n  from: $input_fastq.$dsrc_string.fastq\n  to: $input_fastq.$dsrc_string.fastq.qual\n"
if [ -f $input_fastq.$dsrc_string.fastq.qual ]; then
    printf "$input_fastq.$dsrc_string.fastq.qual already exists (not reproducing it)\n"
else
    $python $xtract_part_fastq_py $input_fastq.$dsrc_string.fastq 3 1> $input_fastq.$dsrc_string.fastq.qual
fi

printf "Compressing quality values with gzip\n  from: $input_fastq.$dsrc_string.fastq.qual\n  to: $input_fastq.$dsrc_string.fastq.qual.gz\n"
if [ -f $input_fastq.$dsrc_string.fastq.qual.gz ]; then
    printf "$input_fastq.$dsrc_string.fastq.qual.gz already exists (not reproducing it)\n"
else
    $gzip -9 -c $input_fastq.$dsrc_string.fastq.qual > $input_fastq.$dsrc_string.fastq.qual.gz
fi

wc -c $input_fastq.$dsrc_string.fastq.qual.gz > $input_fastq.$dsrc_string.fastq.qual.gz.log

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


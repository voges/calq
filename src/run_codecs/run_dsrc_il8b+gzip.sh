#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 input_sam num_threads\n"; exit -1; fi

input_sam=$1
printf "Input SAM file: input_sam\n"
num_threads=$2
printf "Number of threads: $num_threads\n"

if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
dsrc="/project/dna/install/dsrc-2.0/dsrc"
dsrc_string="dsrc-2.0_il8b"
gzip="/usr/bin/gzip"
python="/usr/bin/python"
time="/usr/bin/time"

# Python scripts
sam2fastq_py="/home/voges/git/calq/src/ngstools/sam2fastq.py"
xtract_part_fastq_py="/home/voges/git/calq/src/ngstools/xtract_part_fastq.py"

if [ ! -x $dsrc ]; then printf "Error: Binary file $dsrc is not executable.\n"; exit -1; fi
if [ ! -x $gzip ]; then printf "Error: Binary file $gzip is not executable.\n"; exit -1; fi
if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi
if [ ! -f $sam2fastq_py ]; then printf "Error: Python script $sam2fastq_py is not a regular file.\n"; exit -1; fi
if [ ! -f $xtract_part_fastq_py ]; then printf "Error: Python script $xtract_part_fastq_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "SAM-to-FASTQ conversion\n"
$python $sam2fastq_py $input_sam 1>$input_sam.fastq

printf "Compressing with DSRC\n"
$dsrc c -d3 -q2 -b256 -l -t$num_threads $input_sam.fastq $input_sam.fastq.$dsrc_string

printf "Decompressing with DSRC\n"
$dsrc d -t$num_threads $input_sam.fastq.$dsrc_string $input_sam.fastq.$dsrc_string.fastq

printf "Extracting quality values\n"
$python $xtract_part_fastq_py $input_sam.fastq.$dsrc_string.fastq 3 1> $input_sam.fastq.$dsrc_string.fastq.qual

printf "Compressing quality values with gzip\n"
$gzip -9 -c $input_sam.fastq.$dsrc_string.fastq.qual > $input_sam.fastq.$dsrc_string.fastq.qual.gz

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


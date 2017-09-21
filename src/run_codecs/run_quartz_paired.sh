#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then printf "Usage: $0 input_1_fastq input_2_fastq root\n"; exit -1; fi

input_1_fastq=$1
printf "First input FASTQ file: $input_1_fastq\n"
input_2_fastq=$2
printf "Second input FASTQ file: $input_2_fastq\n"
root=$3
printf "Root: $root\n"

if [ ! -f $input_1_fastq ]; then printf "Error Input FASTQ file $input_1_fastq is not a regular file.\n"; exit -1; fi
if [ ! -f $input_2_fastq ]; then printf "Error Input FASTQ file $input_2_fastq is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
bzip2="/usr/bin/bzip2"
python="/usr/bin/python"
quartz="/project/dna/install/quartz-0.2.2/quartz"
quartz_string="quartz-0.2.2"
quartz_dictionary="/project/dna/resources/quartz_dictionary/dec200.bin.sorted"
time="/usr/bin/time"

# Python scripts
xtract_part_fastq_py="/home/voges/git/calq/src/ngstools/xtract_part_fastq.py"

if [ ! -x $bzip2 ]; then printf "Error: Binary file $bzip2 is not executable.\n"; exit -1; fi
if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -x $quartz ]; then printf "Error: Binary file $quartz is not executable.\n"; exit -1; fi
if [ ! -f $quartz_dictionary ]; then printf "Error: File $quartz_dictionary is not a regular file.\n"; exit -1; fi
if [ ! -x $time ]; then "Error: Binary file $time is not executable.\n"; exit -1; fi
if [ ! -f $xtract_part_fastq_py ]; then printf "Error: Python script xtract_part_fastq_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Running Quartz\n"
cmd="$quartz $quartz_dictionary "$quartz_string" 8 0 $input_1_fastq $input_2_fastq"
$time -v -o $root.$quartz_string.time $cmd &> $root.$quartz_string.log
mv $input_1_fastq.filtered_$quartz_string $input_1_fastq.$quartz_string.fastq
mv $input_2_fastq.filtered_$quartz_string $input_2_fastq.$quartz_string.fastq

printf "Extracting Quartz'd quality values\n"
$python $xtract_part_fastq_py $input_1_fastq.$quartz_string.fastq 3 1> $input_1_fastq.$quartz_string.fastq.qual
$python $xtract_part_fastq_py $input_2_fastq.$quartz_string.fastq 3 1> $input_2_fastq.$quartz_string.fastq.qual

printf "Compressing Quartz'd quality values\n"
$bzip2 $input_1_fastq.$quartz_string.fastq.qual
$bzip2 $input_2_fastq.$quartz_string.fastq.qual

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup ... "
#
printf "Done\n";


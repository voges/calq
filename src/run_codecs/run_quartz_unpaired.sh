#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 1 ]; then printf "Usage: $0 input_fastq\n"; exit -1; fi

input_fastq=$1
printf "Input FASTQ file: $input_fastq\n"

if [ ! -f $input_fastq ]; then printf "Error: Input FASTQ file $input_fastq is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
python="/usr/bin/python"
quartz="/project/dna/install/quartz-0.2.2/quartz"
quartz_string="quartz-0.2.2"
quartz_dictionary="/project/dna/resources/quartz_dictionary/dec200.bin.sorted"
time="/usr/bin/time"

# Python scripts
xtract_part_fastq_py="/home/voges/git/calq/src/ngstools/xtract_part_fastq.py"

if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -x $quartz ]; then printf "Error: Binary file $quartz is not executable.\n"; exit -1; fi
if [ ! -f $quartz_dictionary ]; then printf "Error: File $quartz_dictionary is not a regular file.\n"; exit -1; fi
if [ ! -x $time ]; then "Error: Binary file $time is not executable.\n"; exit -1; fi
if [ ! -f $xtract_part_fastq_py ]; then printf "Error: Python script xtract_part_fastq_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Running Quartz\n"
cmd="$quartz $quartz_dictionary "$quartz_string" 8 0 $input_fastq"
$time -v -o $input_fastq.$quartz_string.time $cmd &> $input_fastq.$quartz_string.log
mv $input_fastq.filtered_$quartz_string $input_fastq.$quartz_string.fastq

printf "Compressing Quartz'd quality values\n"
$python $xtract_part_fastq_py $input_fastq.$quartz_string.fastq 3 1> $input_fastq.$quartz_string.fastq.qual
bzip2 $input_fastq.$quartz_string.fastq.qual
printf "Compressed quality values size: " >> $input_fastq.$quartz_string.log
wc -c $input_fastq.$quartz_string.fastq.qual.bz2 >> $input_fastq.$quartz_string.log
printf "\n" >> $input_fastq.$quartz_string.log

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


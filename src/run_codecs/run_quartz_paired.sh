#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then
    printf "Usage: $0 input_1_fastq input_2_fastq root\n"
    exit -1
fi

input_1_fastq=$1
printf "First input FASTQ file: $input_1_fastq\n"
input_2_fastq=$2
printf "Second input FASTQ file: $input_2_fastq\n"
root=$3
printf "Root: $root\n"

printf "Checking first input FASTQ file $input_1_fastq ... "
if [ ! -f $input_1_fastq ]; then printf "did not find first input FASTQ file: $input_1_fastq\n"; exit -1; fi
printf "OK\n"

printf "Checking second input FASTQ file $input_2_fastq ... "
if [ ! -f $input_2_fastq ]; then printf "did not find second input FASTQ file: $input_2_fastq\n"; exit -1; fi
printf "OK\n"

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

printf "Checking executables ... "
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -x $quartz ]; then printf "did not find $quartz\n"; exit -1; fi
if [ ! -e $quartz_dictionary ]; then printf "did not find $quartz_dictionary\n"; exit -1; fi
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -e $xtract_part_fastq_py ]; then printf "did not find $xtract_part_fastq_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Running Quartz ... "
cmd="$quartz $quartz_dictionary "$quartz_string" 8 0 $input_1_fastq $input_2_fastq"
$time -v -o $root.$quartz_string.time $cmd &> $root.$quartz_string.log &
mv $input_1_fastq.filtered_$quartz_string $input_1_fastq.$quartz_string.fastq
mv $input_2_fastq.filtered_$quartz_string $input_2_fastq.$quartz_string.fastq
printf "OK\n"

printf "Compressing Quartz'd quality values ... "
$python $xtract_part_fastq_py $input_1_fastq.$quartz_string.fastq 3 1> $input_1_fastq.$quartz_string.fastq.qual
$python $xtract_part_fastq_py $input_2_fastq.$quartz_string.fastq 3 1> $input_2_fastq.$quartz_string.fastq.qual
bzip2 $input_1_fastq.$quartz_string.fastq.qual
bzip2 $input_2_fastq.$quartz_string.fastq.qual
printf "Compressed quality values size: " >> $input_1_fastq.$quartz_string.log
printf "Compressed quality values size: " >> $input_2_fastq.$quartz_string.log
wc -c $input_1_fastq.$quartz_string.fastq.qual.bz2 >> $input_1_fastq.$quartz_string.log
wc -c $input_2_fastq.$quartz_string.fastq.qual.bz2 >> $input_2_fastq.$quartz_string.log
printf "\n" >> $input_1_fastq.$quartz_string.log
printf "\n" >> $input_2_fastq.$quartz_string.log
printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


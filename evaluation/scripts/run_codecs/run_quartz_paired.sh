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

pgrep="/usr/bin/pgrep"
python="/usr/bin/python"
time="/usr/bin/time"

ps_mem_py="/home/voges/git/calq/evaluation/scripts/ps_mem/ps_mem.py"
xtract_qual_fastq_py="/home/voges/git/ngstools/xtract_qual_fastq.py"

quartz="/project/dna/install/quartz-0.2.2/quartz"
quartz_string="quartz-0.2.2"
quartz_dictionary="/project/dna/resources/quartz_dictionary/dec200.bin.sorted"

printf "Checking executables ... "
if [ ! -x $pgrep ]; then printf "did not find $pgrep\n"; exit -1; fi
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -e $ps_mem_py ]; then printf "did not find $ps_mem_py\n"; exit -1; fi
if [ ! -e $xtract_qual_fastq_py ]; then printf "did not find $xtract_qual_fastq_py\n"; exit -1; fi
if [ ! -x $quartz ]; then printf "did not find $quartz\n"; exit -1; fi
if [ ! -e $quartz_dictionary ]; then printf "did not find $quartz_dictionary\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Running Quartz ... "
cmd="$quartz $quartz_dictionary "$quartz_string" 8 0 $input_1_fastq $input_2_fastq"
$time -v -o $root.$quartz_string.time $cmd &> $root.$quartz_string.log &
time_pid=$!
cmd_pid=$($pgrep -P $time_pid)
printf "Command being traced: \"$cmd\"\n" > $root.$quartz_string.mem
$python $ps_mem_py -t -w 1 --swap -p $cmd_pid >> $root.$quartz_string.mem
mv $input_1_fastq.filtered_$quartz_string $input_1_fastq.$quartz_string.fastq
mv $input_2_fastq.filtered_$quartz_string $input_2_fastq.$quartz_string.fastq
printf "OK\n"

printf "Compressing Quartz'd quality values ... "
$python $xtract_qual_fastq_py $input_1_fastq.$quartz_string.fastq 2> $input_1_fastq.$quartz_string.fastq.qual
$python $xtract_qual_fastq_py $input_2_fastq.$quartz_string.fastq 2> $input_2_fastq.$quartz_string.fastq.qual
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


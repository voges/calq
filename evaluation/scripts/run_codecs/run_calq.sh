#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then
    printf "Usage: $0 input_sam polyploidy block_size\n"
    exit -1
fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
polyploidy=$2
printf "Polyploidy: $polyploidy\n"
block_size=$3
printf "Block size: $block_size\n"

printf "Checking input SAM file $input_sam ... "
if [ ! -f $input_sam ]; then printf "did not find input SAM file: $input_sam\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

pgrep="/usr/bin/pgrep"
python="/usr/bin/python"
time="/usr/bin/time"

replace_qual_sam_py="/home/voges/git/ngstools/replace_qual_sam.py"
ps_mem_py="/home/voges/git/calq/evaluation/scripts/ps_mem/ps_mem.py"

calq="/project/dna/install/calq-d3965e2/calq"
calq_string="calq-d3965e2"

printf "Checking executables ... "
if [ ! -x $pgrep ]; then printf "did not find $pgrep\n"; exit -1; fi
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -e $replace_qual_sam_py ]; then printf "did not find $replace_qual_sam_py\n"; exit -1; fi
if [ ! -e $ps_mem_py ]; then printf "did not find $ps_mem_py\n"; exit -1; fi
if [ ! -x $calq ]; then printf "did not find $calq\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                          Compress and decompress                            #
###############################################################################

printf "Running CALQ encoder ... "
cmd="$calq -f -q Illumina-1.8+ -p $polyploidy -b $block_size $input_sam -o $input_sam.cq"
$time -v -o $input_sam.$calq_string.enc.time $cmd &> $input_sam.$calq_string.enc.log & time_pid=$!
cmd_pid=$($pgrep -P $time_pid)
printf "Command being traced: \"$cmd\"\n" > $input_sam.$calq_string.enc.mem
$python $ps_mem_py -t -w 1 --swap -p $cmd_pid >> $input_sam.$calq_string.enc.mem
mv $input_sam.cq $input_sam.$calq_string
printf "OK\n"

printf "Running CALQ decoder ... "
cmd="$calq -f -d -s $input_sam $input_sam.$calq_string -o $input_sam.$calq_string.qual"
$time -v -o $input_sam.$calq_string.dec.time $cmd &> $input_sam.$calq_string.dec.log & time_pid=$!
cmd_pid=$($pgrep -P $time_pid)
printf "Command being traced: \"$cmd\"\n" > $input_sam.$calq_string.dec.mem
$python $ps_mem_py -t -w 1 --swap -p $cmd_pid >> $input_sam.$calq_string.dec.mem
printf "OK\n"

#printf "Constructing SAM file with reconstruced quality values ... "
#$python $replace_qual_sam_py $input_sam $input_sam.$calq_string.qual &> $input_sam.run_calq.log
#mv $input_sam.new_qual.sam $input_sam.$calq_string.sam
#printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
rm $input_sam.$calq_string.qual
printf "OK\n";


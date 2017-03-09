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

time="/usr/bin/time"
python="/usr/bin/python"
calq="/project/dna/install/calq-d3965e2/calq"
calq_string="calq-d3965e2"
replace_qual_sam_py="/home/voges/git/ngstools/replace_qual_sam.py"

printf "Checking executables ... "
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -x $calq ]; then printf "did not find $calq\n"; exit -1; fi
if [ ! -e $replace_qual_sam_py ]; then printf "did not find $replace_qual_sam_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                          Compress and decompress                            #
###############################################################################

printf "[1/3] Running CALQ encoder ... "
$time -v -o $input_sam.$calq_string.enc.time $calq -f -q Illumina-1.8+ -p $polyploidy -b $block_size $input_sam -o $input_sam.cq &> $input_sam.$calq_string.enc.log
printf "OK\n"

printf "[2/3] Running CALQ decoder ... "
$time -v -o $input_sam.$calq_string.dec.time $calq -f -d -s $input_sam $input_sam.cq &> $input_sam.$calq_string.dec.log
printf "OK\n"

printf "[3/3] Constructing SAM file with reconstruced quality values ... "
mv $input_sam.cq $input_sam.$calq_string
mv $input_sam.cq.qual $input_sam.$calq_string.qual
$python $replace_qual_sam_py $input_sam $input_sam.$calq_string.qual &> $input_sam.run_calq.log
mv $input_sam.new_qual.sam $input_sam.$calq_string.sam
printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


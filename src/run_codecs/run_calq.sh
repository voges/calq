#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 4 ]; then printf "Usage: $0 input_sam qual_type polyploidy block_size\n"; exit -1; fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
qual_type=$2
printf "Quality value type: $qual_type\n"
polyploidy=$3
printf "Polyploidy: $polyploidy\n"
block_size=$4
printf "Block size: $block_size\n"

if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
calq="/project/dna/install/calq-616e792/calq"
calq_string="calq-616e792"
python="/usr/bin/python"
time="/usr/bin/time"

# Python scripts
replace_qual_sam_py="/home/voges/git/calq/src/ngstools/replace_qual_sam.py"

if [ ! -x $calq ]; then printf "Error: Binary file $calq is not executable.\n"; exit -1; fi
if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi
if [ ! -f $replace_qual_sam_py ]; then printf "Error: Python script $replace_qual_sam_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                          Compress and decompress                            #
###############################################################################

printf "Compressing with CALQ\n  from: $input_sam\n  to: $input_sam.$calq_string\n"
cmd="$calq -q $qual_type -p $polyploidy -b $block_size $input_sam -o $input_sam.$calq_string"
if [ -f $input_sam.$calq_string ]; then
    printf "$input_sam.$calq_string already exists (not reproducing it)\n"
else
    $time -v -o $input_sam.$calq_string.enc.time $cmd &> $input_sam.$calq_string.enc.log
fi

printf "Decompressing with CALQ\n  from: $input_sam.$calq_string\n  to: $input_sam.$calq_string.qual\n"
cmd="$calq -d -s $input_sam $input_sam.$calq_string -o $input_sam.$calq_string.qual"
if [ -f $input_sam.$calq_string.qual ]; then
    printf "$input_sam.$calq_string.qual already exists (not reproducing it)\n"
else
    $time -v -o $input_sam.$calq_string.dec.time $cmd &> $input_sam.$calq_string.dec.log
fi

printf "Constructing SAM file with reconstruced quality values\n  from: $input_sam.$calq_string.qual\n  to: $input_sam.$calq_string.sam"
if [ -f $input_sam.$calq_string.sam ]; then
    printf "$input_sam.$calq_string.sam already exists (not reproducing it)\n"
else
    $python $replace_qual_sam_py $input_sam $input_sam.$calq_string.qual 1> $input_sam.$calq_string.sam
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


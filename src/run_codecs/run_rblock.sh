#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 input_sam r (with r={5,10,20,40,80,160,320,640})\n"; exit -1; fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
r=$1
printf "r: $r\n" # r={5,10,20,40,80,160,320,640}

if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

prblock_compress="/project/dna/prog/libCSAM-da36a12/CompressQual"
prblock_decompress="/project/dna/prog/libCSAM-da36a12/DecompressQual"
rblock_string="rblock"
time="/usr/bin/time"

if [ ! -x $prblock_compress ]; then printf "Error: Binary file $prblock_compress is not executable.\n"; exit -1; fi
if [ ! -x $prblock_decompress ]; then printf "Error: Binary file $prblock_decompress is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Compressing quality values with r=$r\n"
cmd="$prblock_compress $input_sam -q 2 -m 1 -l $r"
$time -v -o $input_sam.$rblock_string.enc.time $cmd &> $input_sam.$rblock_string.enc.log

printf "Decompressing quality values\n"
cmd="$prblock_decompress $input_sam.cqual"
$time -v -o $input_sam.$rblock_string.enc.time $cmd &> $input_sam.$rblock_string.enc.log

mv $input_sam.cqual $input_sam.qual.$rblock_string
mv $input_sam.cqual.qual $input_sam.qual.$rblock_string.qual

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 input_sam p (with p={1,2,4,8,16,32})\n"; exit -1; fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
p=$2
printf "p: $p\n" # p={1,2,4,8,16,32}

if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

prblock_compress="/project/dna/prog/libCSAM-da36a12/CompressQual"
prblock_decompress="/project/dna/prog/libCSAM-da36a12/DecompressQual"
pblock_string="pblock"
time="/usr/bin/time"

if [ ! -x $prblock_compress ]; then printf "Error: Binary file $prblock_compress is not executable.\n"; exit -1; fi
if [ ! -x $prblock_decompress ]; then printf "Error: Binary file $prblock_decompress is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Compressing quality values with p=$p\n"
$prblock_compress $input_sam -q 1 -m 1 -l $p

printf "Decompressing quality values\n"
$prblock_decompress $input_sam.cqual

mv $input_sam.cqual $input_sam.qual.$pblock_string
mv $input_sam.qual $input_sam.qual.$pblock_string.qual
wc -c $input_sam.qual.$pblock_string > $input_sam.qual.$pblock_string.log

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


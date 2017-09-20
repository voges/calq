#! /bin/bash

if [ "$#" -ne 1 ]; then printf "Usage: $0 file.sam\n"; exit -1; fi
input_sam=$1
printf "SAM file: $input_sam\n"
if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi
samtools="/project/dna/install/samtools-1.3/bin/samtools"
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
$samtools depth $1 | awk '{sum+=$3} END { print "Coverage: ",sum/NR }'

#!/bin/bash

if [ "$#" -ne 1 ]; then printf "Usage: $0 file.[fastq|fq] 1>file.[fasta|fa]\n"; exit -1; fi
input_fastq=$1
if [ ! -f $input_fastq ]; then printf "Error: Input FASTQ file $input_fastq is not a regular file.\n"; exit -1; fi
cat $input_fastq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}'

#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 file.[fa|fasta]"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .fa or .fasta

install_path="/project/dna/install"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"

if [ -f $root.dict ]; then
    date; echo "FASTA dict file $root.dict already exists. Not reproducing it."
else
    date; echo "Constructing FASTA dict file $root.dict"
    java -jar $picard_jar CreateSequenceDictionary R=$1 O=$root.dict
fi


#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 file.[fa|fasta]"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .fa or .fasta

install_path="/project/dna/install"
samtools="$install_path/samtools-1.3/bin/samtools"

if [ -f $1.fai ]; then
    echo "FASTA index file $1.fai already exists. Not reproducing it."
else
    echo "Constructing FASTA index file $1.fai"
    $samtools faidx $1
fi


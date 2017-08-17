#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file.[fa|fasta] chromosome"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .fa or .fasta
chromosome=$2

install_path="/project/dna/install"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
samtools="$install_path/samtools-1.3/bin/samtools"

echo "Creating FASTA index file $1.fai"
if [ -f $1.fai ]; then
    echo "FASTA index file $1.fai already exists. Not reproducing it."
else
    $samtools faidx $1
fi

echo "Writing chromosome $chromosome from FASTA file $1 to FASTA file $root.$chromosome.fasta"
if [ -f $root.$chromosome.fasta ]; then
    echo "FASTA file $root.$chromosome.fasta already exists. Not reproducing it."
else
    $samtools faidx $1 $chromosome 1> $root.$chromosome.fasta
fi

echo "Creating FASTA index file $root.$chromosome.fasta.fai"
if [ -f $root.$chromosome.fasta.fai ]; then
    echo "FASTA index file $root.$chromosome.fasta.fai already exists. Not reproducing it."
else
    $samtools faidx $root.$chromosome.fasta
fi

echo "Creating FASTA dict file $root.$chromosome.dict"
if [ -f $root.$chromosome.dict ]; then
    echo "FASTA dict file $root.$chromosome.dict already exists. Not reproducing it."
else
    java -jar $picard_jar CreateSequenceDictionary R=$root.$chromosome.fasta O=$root.$chromosome.dict
fi


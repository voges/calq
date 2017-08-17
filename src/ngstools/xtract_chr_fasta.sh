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

if [ -f $1.fai ]; then
    echo "FASTA index file $1.fai already exists. Not reproducing it."
    #touch $1.fai
else
    echo "Constructing FASTA index file $1.fai"
    $samtools faidx $1
fi

echo "Extracting chromosome $chromosome from FASTA file $1"
$samtools faidx $1 $chromosome 1> $root.$chromosome.fasta
echo "Wrote chromosome $chromosome to FASTA file $root.$chromosome.fasta"
echo "Constructing FASTA index file $root.$chromosome.fasta.fai"
$samtools faidx $root.$chromosome.fasta
echo "Finished constructing FASTA index file $root.$chromosome.fasta.fai"
echo "Constructing FASTA dict file $root.$chromosome.dict"
java -jar $picard_jar CreateSequenceDictionary R=$root.$chromosome.fasta O=$root.$chromosome.dict
echo "Finished constructing FASTA index file $root.$chromosome.fasta.fai"


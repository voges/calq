#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file.[sam|bam] chromosome"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .sam/.bam
chromosome=$2

install_path="/project/dna/install"
samtools="$install_path/samtools-1.3/bin/samtools"

if [[ $1 == *.sam ]]; then
    echo "Converting SAM file $1 to BAM file $root.bam"
    if [ -f $root.bam ]; then
        echo "BAM file $root.bam already exists. Not reproducing it."
    else
        $samtools view -bh $1 > $root.bam
    fi
fi

echo "Constructing BAM index file $root.bai"
if [ -f $root.bai ]; then
    echo "BAM index file $root.bai already exists. Not reproducing it."
else
    $samtools index $root.bam $root.bai
fi

echo "Writing chromosome $chromosome from $root.bam to $root.$chromosome.sam"
if [ -f $root.$chromosome.sam ]; then
    echo "SAM file $root.$chromosome.sam already exists. Not reproducing it."
else
    $samtools view -h $root.bam $chromosome > $root.$chromosome.sam
fi


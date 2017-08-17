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
    $samtools view -bh $1 > $root.bam
fi

if [ -f $root.bai ]; then
    echo "BAM index file $root.bai already exists. Not reproducing it."
    touch $root.bai
else
    echo "Constructing BAM index file $root.bai"
    $samtools index $root.bam $root.bai
fi

echo "Extracting chromosome $chromosome from $root.bam"
$samtools view -h $root.bam $chromosome > $root.$chromosome.sam
echo "Wrote chromosome $chromosome to $root.$chromosome.sam"


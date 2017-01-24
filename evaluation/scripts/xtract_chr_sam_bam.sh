#!/bin/bash

###############################################################################
#               Extract chromosome from SAM or BAM file                       #
###############################################################################

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file.[sam|bam] chromosome"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .sam/.bam
chromosome=$2

install_path="/project/dna/install"
samtools="$install_path/samtools-1.3/bin/samtools"

if [[ $1 == *.sam ]]; then
    date; echo "Converting SAM file $1 to BAM file $root.bam"
    $samtools view -bh $1 > $root.bam
fi

if [ -f $root.bai ]; then
    date; echo "BAM index file $root.bai already exists. Not reproducing it."
    touch $root.bai
else
    date; echo "Constructing BAM index file $root.bai"
    $samtools index $root.bam $root.bai
fi

date; echo "Extracting chromosome $chromosome from $root.bam"
$samtools view -h $root.bam $chromosome > $root.$chromosome.sam
date; echo "Wrote chromosome $chromosome to $root.$chromosome.sam"

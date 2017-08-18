#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 file.[sam|bam] chromosome\n"; exit -1; fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .sam/.bam
chromosome=$2

printf "Input file: $1\n"
printf "Chromosome: $chromosome\n"

if [ ! -f $1 ]; then printf "Error: Input file $1 is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

install_path="/project/dna/install"
samtools="$install_path/samtools-1.3/bin/samtools"

if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi

###############################################################################
#                                 Extraction                                  #
###############################################################################

if [[ $1 == *.sam ]]; then
    printf "SAM-to-BAM conversion\n  from: $1\n  to: $root.bam\n"
    if [ -f $root.bam ]; then
        printf "$root.bam already exists (not reproducing it)\n"
    else
        $samtools view -bh $1 1> $root.bam
    fi
fi

printf "Constructing BAM index file: $root.bai\n"
if [ -f $root.bai ]; then
    printf "$root.bai already exists (not reproducing it)\n"
else
    $samtools index $root.bam $root.bai
fi

printf "Extracting\n  chromosome: $chromosome\n  from: $root.bam\n  to: $root.$chromosome.sam\n"
if [ -f $root.$chromosome.sam ]; then
    printf "$root.$chromosome.sam already exists (not reproducing it)\n"
else
    $samtools view -h $root.bam $chromosome 1> $root.$chromosome.sam
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


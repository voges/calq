#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 file.[fa|fasta] chromosome\n"; exit -1; fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .fa or .fasta
chromosome=$2

printf "Input FASTA file: $1\n"
printf "Chromosome: $chromosome\n"

if [ ! -f $1 ]; then printf "Error: Input FASTA file $1 is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

install_path="/project/dna/install"
java="/usr/bin/java"
java_opts=""
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
samtools="$install_path/samtools-1.3/bin/samtools"

if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $picard_jar ]; then printf "Error: JAR file $picard_jar is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi

###############################################################################
#                                 Extraction                                  #
###############################################################################

printf "Creating FASTA index file: $1.fai\n"
if [ -f $1.fai ]; then
    printf "$1.fai already exists (not reproducing it)\n"
else
    $samtools faidx $1
fi

printf "Extracting\n  chromosome: $chromosome\n  from: $1\n  to: $root.$chromosome.fasta\n"
if [ -f $root.$chromosome.fasta ]; then
    printf "$root.$chromosome.fasta already exists (not reproducing it)\n"
else
    $samtools faidx $1 $chromosome 1> $root.$chromosome.fasta
fi

printf "Creating FASTA index file: $root.$chromosome.fasta.fai\n"
if [ -f $root.$chromosome.fasta.fai ]; then
    printf "$root.$chromosome.fasta.fai already exists (not reproducing it)\n"
else
    $samtools faidx $root.$chromosome.fasta
fi

printf "Creating FASTA dict file: $root.$chromosome.dict\n"
if [ -f $root.$chromosome.dict ]; then
    printf "$root.$chromosome.dict already exists (not reproducing it)\n"
else
    $java $java_opts -jar $picard_jar CreateSequenceDictionary \
        R=$root.$chromosome.fasta \
        O=$root.$chromosome.dict
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


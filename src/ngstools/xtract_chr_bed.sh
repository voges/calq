#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 file.bed chromosome\n"; exit -1; fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .bed
chromosome=$2

printf "Input BED file: $1\n"
printf "Chromosome: $chromosome\n"

if [ ! -f $1 ]; then printf "Error: Input BED file $1 is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

grep="/usr/bin/grep"

if [ ! -x $grep ]; then printf "Error: Binary file $grep is not executable.\n"; exit -1; fi

###############################################################################
#                                 Extraction                                  #
###############################################################################

printf "Extracting\n  chromosome: $chromosome\n  from: $1\n  to: $root.$chromosome.bed\n"
if [ -f $root.$chromosome.bed ]; then
    printf "$root.$chromosome.bed already exists (not reproducing it)\n"
else
    grep -P "^$chromosome" $1 1> $root.$chromosome.bed
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


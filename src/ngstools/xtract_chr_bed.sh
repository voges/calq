#!/bin/bash

###############################################################################
#                   Extract chromosome from BED file                          #
###############################################################################

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file.bed chromosome"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .bed
chromosome=$2

if [ -f $root.$chromosome.bed ]; then
    date; echo "BED file $root.$chromosome.bed already exists. Aborting."
else
    date; echo "Extracting chromosome $chromosome from BED file $1"
    grep -P "^$chromosome" $1 > $root.$chromosome.bed
    date; echo "Wrote chromosome $chromosome data to BED file $root.$chromosome.bed"
fi

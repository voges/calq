#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file.bed chromosome"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .bed
chromosome=$2

echo "Writing chromosome $chromosome from BED file $1 to BED file $root.$chromosome.bed"
if [ -f $root.$chromosome.bed ]; then
    echo "BED file $root.$chromosome.bed already exists. Not reproducing it."
else
    grep -P "^$chromosome" $1 > $root.$chromosome.bed
fi


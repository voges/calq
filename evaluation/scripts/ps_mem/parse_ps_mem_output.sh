#!/bin/bash

# Parse ps_mem.py output produced by e.g. the following command:
#   $ ps_mem.py -p 17449 -t -w 1 --swap
# The output of ps_mem.py would look similar to the following snippet:
#   RAM: 119123  Swap: 754370
#   RAM: 11951616  Swap: 23423
#   RAM: 3121951616  Swap: 9675450

if [ "$#" -ne 1 ]; then
    printf "Usage: $0 input_mem\n"
    exit -1
fi

# Delete the first and last two lines of the input file
head -n -2 $1 | tail -n +3 > $1.tmp

printf "RAM usage:\n----------\n"
printf "average: "
awk '{ sum += $2; n++ } END { if (n > 0) print sum/n; }' $1.tmp
printf "min: "
awk 'BEGIN {min = 1000000000000} {if ($2<min) min=$2} END {print min}' $1.tmp
printf "max: "
awk 'BEGIN {max = 0} {if ($2>max) max=$2} END {print max}' $1.tmp

printf "\nSwap usage:\n----------\n"
printf "average: "
awk '{sum += $4; n++} END { if (n > 0) print sum/n; }' $1.tmp
printf "min: "
awk 'BEGIN {min = 1000000000000} {if ($4<min) min=$4} END {print min}' $1.tmp
printf "max: "
awk 'BEGIN {max = 0} {if ($4>max) max=$4} END {print max}' $1.tmp

rm -f $1.tmp


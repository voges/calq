#! /bin/bash
samtools depth $1 | awk '{sum+=$3} END { print "Coverage: ",sum/NR }'

#! /bin/bash
/project/dna/bin/samtools depth $1 | awk '{sum+=$3} END { print "Coverage: ",sum/NR }'

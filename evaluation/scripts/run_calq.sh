#!/bin/bash
set -x
/project/dna/bin/calq-d3965e2 -q Illumina-1.8+ -p 1 -b 10000 -f $1 -o $1.cq 2>&1 | tee $1.calq-d3965e2.enc.log
/project/dna/bin/calq-d3965e2 -d -s $1 -f $1.cq 2>&1 | tee $1.calq-d3965e2.dec.log
mv $1.cq $1.calq-d3965e2
mv $1.cq.qual $1.calq-d3965e2.qual
/home/voges/git/calq/evaluation/scripts/replace_qual_sam.py $1 $1.calq-d3965e2.qual
mv $1.new_qual.sam $1.calq-d3965e2.sam

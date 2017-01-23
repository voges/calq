#!/bin/bash
set -x
/project/dna/prog/quartz-0.2.2/quartz /project/dna/resources/quartz_dictionary/dec200.bin.sorted "quartz-0.2.2" 8 0 $1_1.fastq $1_2.fastq 2>&1 | tee $1.quartz-0.2.2.enc.log
mv $1.filtered_quartz-0.2.2 $1.quartz-0.2.2.fastq
/home/voges/git/calq/evaluation/scripts/xtract_qual_fastq.py $1.quartz-0.2.2.fastq 2> $1.quartz-0.2.2.fastq.qual
bzip2 $1.quartz-0.2.2.fastq.qual
wc -c $1.quartz-0.2.2.fastq.qual.bz2

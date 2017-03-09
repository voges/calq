#!/bin/bash


if [ "$#" -ne 1 ]; then
    echo "Usage: $0 fastq_file"
    exit -1
fi

date
set -x

fastq_file=$1
install_path="/project/dna/install"
quartz="$install_path/quartz-0.2.2/quartz"
quartz_string="quartz-0.2.2"
quartz_dictionary="/project/dna/resources/quartz_dictionary/dec200.bin.sorted"

$quartz $quartz_dictionary "$quartz_string" 8 0 $fastq_file 2>&1 | tee $fastq_file.$quartz_string.log
mv $fastq_file.filtered_$quartz_string $fastq_file.$quartz_string.fastq
/home/voges/git/calq/evaluation/scripts/xtract_qual_fastq.py $fastq_file.$quartz_string.fastq 2> $fastq_file.$quartz_string.fastq.qual
bzip2 $fastq_file.$quartz_string.fastq.qual
wc -c $fastq_file.$quartz_string.fastq.qual.bz2

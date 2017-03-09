#!/bin/bash


if [ "$#" -ne 1 ]; then
    echo "Usage: $0 root"
    exit -1
fi

date
set -x

root=$1
install_path="/project/dna/install"
quartz="$install_path/quartz-0.2.2/quartz"
quartz_string="quartz-0.2.2"
quartz_dictionary="/project/dna/resources/quartz_dictionary/dec200.bin.sorted"

$quartz $quartz_dictionary "$quartz_string" 8 0 $root\_1.fastq $root\_2.fastq 2>&1 | tee $root.$quartz_string.log
mv $root\_1.fastq.filtered_$quartz_string $root\_1.fastq.$quartz_string.fastq
mv $root\_2.fastq.filtered_$quartz_string $root\_2.fastq.$quartz_string.fastq
/home/voges/git/calq/evaluation/scripts/xtract_qual_fastq.py $root\_1.fastq.$quartz_string.fastq 2> $root\_1.fastq.$quartz_string.fastq.qual
/home/voges/git/calq/evaluation/scripts/xtract_qual_fastq.py $root\_2.fastq.$quartz_string.fastq 2> $root\_2.fastq.$quartz_string.fastq.qual
bzip2 $root\_1.fastq.$quartz_string.fastq.qual
bzip2 $root\_2.fastq.$quartz_string.fastq.qual
wc -c $root\_1.fastq.$quartz_string.fastq.qual.bz2
wc -c $root\_2.fastq.$quartz_string.fastq.qual.bz2

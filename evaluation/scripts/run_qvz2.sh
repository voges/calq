#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 fastq_file T"
    exit -1
fi

date
set -x

fastq_file=$1
T=$2
install_path="/project/dna/install"
qvz2="$install_path/qvz2-d5383c6/qvz"
qvz2_string="qvz2-d5383c6"

$qvz2 -t $T -v -u $fastq_file.$qvz2_string-t$T.qual $fastq_file.qual $fastq_file.$qvz2_string-t$T 2>&1 | tee $fastq_file.$qvz2_string-t$T.log
/home/voges/git/calq/evaluation/scripts/replace_qual_sam.py $fastq_file $fastq_file.$qvz2_string-t$T.qual
mv $fastq_file.new_qual.sam $fastq_file.$qvz2_string-t$T.sam

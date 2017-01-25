#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 sam_file polyploidy"
    exit -1
fi

date
set -x

sam_file=$1
polyploidy=$2
install_path="/project/dna/install"
calq="$install_path/calq-d3965e2/calq"
calq_string="calq-d3965e2"

$calq -q Illumina-1.8+ -p $polyploidy -b 10000 -f $sam_file -o $sam_file.cq 2>&1 | tee $sam_file.$calq_string.enc.log
$calq -d -s $sam_file -f $sam_file.cq 2>&1 | tee $sam_file.$calq_string.dec.log
mv $sam_file.cq $sam_file.$calq_string
mv $sam_file.cq.qual $sam_file.$calq_string.qual
/home/voges/git/calq/evaluation/scripts/replace_qual_sam.py $sam_file $sam_file.$calq_string.qual
mv $sam_file.new_qual.sam $sam_file.$calq_string.sam

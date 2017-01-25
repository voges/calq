#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 bam_file ref_file"
    exit -1
fi

date
set -x

bam_file=$1
ref_file=$2
install_path="/project/dna/install"
crumble="$install_path/crumble-0.5/crumble"
crumble_string="crumble-0.5"
scramble="$install_path/io_lib-1.14.6/bin/scramble"
scramble_string="scramble-1.14.6"
cram_size="$install_path/io_lib-1.14.6/bin/cram_size"

$crumble -v -1 $bam_file $bam_file.$crumble_string-1.bam 2>&1 | tee $bam_file.$crumble_string-1.log
$scramble -r $ref_file -t16 $bam_file.$crumble_string-1.bam $bam_file.$crumble_string-1.bam.cram 2>&1 | tee $bam_file.$crumble_string-1.bam.$scramble_string-1.log
$cram_size $bam_file.$crumble_string-1.bam.cram 2>&1 | tee $bam_file.$crumble_string-1.bam.$scramble_string.cram_size
mv $bam_file.$crumble_string-1.bam.cram $bam_file.$crumble_string-1.bam.$scramble_string

$crumble -v -9 $bam_file $bam_file.$crumble_string-9.bam 2>&1 | tee $bam_file.$crumble_string-9.log
$scramble -r $ref_file -t16 $bam_file.$crumble_string-9.bam $bam_file.$crumble_string-9.bam.cram 2>&1 | tee $bam_file.$crumble_string-9.bam.$scramble_string.log
$cram_size $bam_file.$crumble_string-1.bam.cram 2>&1 | tee $bam_file.$crumble_string-1.bam.$scramble_string.cram_size
mv $bam_file.$crumble_string-9.bam.cram $bam_file.$crumble_string-9.bam.$scramble_string

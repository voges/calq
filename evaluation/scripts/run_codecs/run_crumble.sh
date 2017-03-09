#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then
    printf "Usage: $0 input_bam ref_fasta\n"
    exit -1
fi

input_bam=$1
printf "Input BAM file: $input_bam\n"
ref_fasta=$2
printf "Reference FASTA file: $ref_fasta\n"

printf "Checking input BAM file $input_bam ... "
if [ ! -f $input_bam ]; then printf "did not find input BAM file: $input_bam\n"; exit -1; fi
printf "OK\n"

printf "Checking reference FASTA file $ref_fasta ... "
if [ ! -f $ref_fasta ]; then printf "did not find reference FASTA file: $ref_fasta\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

time="/usr/bin/time"
crumble="/project/dna/install/crumble-0.5/crumble"
crumble_string="crumble-0.5"
scramble="/project/dna/install/io_lib-1.14.6/bin/scramble"
scramble_string="scramble-1.14.6"
cram_size="/project/dna/install/io_lib-1.14.6/bin/cram_size"

printf "Checking executables ... "
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -x $crumble ]; then printf "did not find $crumble\n"; exit -1; fi
if [ ! -x $scramble ]; then printf "did not find $scramble\n"; exit -1; fi
if [ ! -x $cram_size ]; then printf "did not find $cram_size\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

$time -v -o $input_bam.$crumble_string-1.time $crumble -v -1 $input_bam $input_bam.$crumble_string-1.bam &> $input_bam.$crumble_string-1.log
$time -v -o $input_bam.$crumble_string-1.bam.$scramble_string-1.time $scramble -r $ref_fasta -t16 $input_bam.$crumble_string-1.bam $input_bam.$crumble_string-1.bam.cram &> $input_bam.$crumble_string-1.bam.$scramble_string-1.log
$cram_size $input_bam.$crumble_string-1.bam.cram &> $input_bam.$crumble_string-1.bam.$scramble_string.cram_size
mv $input_bam.$crumble_string-1.bam.cram $input_bam.$crumble_string-1.bam.$scramble_string

$time -v -o $input_bam.$crumble_string-9.time $crumble -v -9 $input_bam $input_bam.$crumble_string-9.bam &> $input_bam.$crumble_string-9.log
$time -v -o $input_bam.$crumble_string-9.bam.$scramble_string-1.time $scramble -r $ref_fasta -t16 $input_bam.$crumble_string-9.bam $input_bam.$crumble_string-9.bam.cram &> $input_bam.$crumble_string-9.bam.$scramble_string.log
$cram_size $input_bam.$crumble_string-9.bam.cram &> $input_bam.$crumble_string-9.bam.$scramble_string.cram_size
mv $input_bam.$crumble_string-9.bam.cram $input_bam.$crumble_string-9.bam.$scramble_string

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


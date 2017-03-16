#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then
    printf "Usage: $0 input_bam ref_fasta num_threads\n"
    exit -1
fi

input_bam=$1
printf "Input BAM file: $input_bam\n"
ref_fasta=$2
printf "Reference FASTA file: $ref_fasta\n"
num_threads=$3
printf "Number of threads: $num_threads\n"

printf "Checking input BAM file $input_bam ... "
if [ ! -f $input_bam ]; then printf "did not find input BAM file: $input_bam\n"; exit -1; fi
printf "OK\n"

printf "Checking reference FASTA file $ref_fasta ... "
if [ ! -f $ref_fasta ]; then printf "did not find reference FASTA file: $ref_fasta\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

pgrep="/usr/bin/pgrep"
python="/usr/bin/python"
time="/usr/bin/time"
crumble="/project/dna/install/crumble-0.5/crumble"
crumble_string="crumble-0.5"
scramble="/project/dna/install/io_lib-1.14.6/bin/scramble"
scramble_string="scramble-1.14.6"
cram_size="/project/dna/install/io_lib-1.14.6/bin/cram_size"
samtools="/project/dna/install/samtools-1.3/bin/samtools"

printf "Checking executables ... "
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -x $crumble ]; then printf "did not find $crumble\n"; exit -1; fi
if [ ! -x $scramble ]; then printf "did not find $scramble\n"; exit -1; fi
if [ ! -x $cram_size ]; then printf "did not find $cram_size\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "did not find $samtools\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Checking if FASTA index file exists ... "
if [ -f $ref_fasta.fai ]; then
    printf "FASTA index file $ref_fasta.fai already exists. Not reproducing it.\n"
else
    printf "FASTA index file $ref_fasta.fai does not exist. Constructing FASTA index file $ref_fasta.fai ... "
    $samtools faidx $ref_fasta
    printf "OK\n"
fi

printf "Running Crumble BAM-to-BAM encoder with compression level 1 ... "
cmd="$crumble -v -1 $input_bam $input_bam.$crumble_string-1.bam"
$time -v -o $input_bam.$crumble_string-1.time $cmd &> $input_bam.$crumble_string-1.log & time_pid=$!
cmd_pid=$($pgrep -P $time_pid)
printf "Command being traced: \"$cmd\"\n" > $input_bam.$crumble_string-1.mem
$python $ps_mem_py -t -w 1 --swap -p $cmd_pid >> $input_bam.$crumble_string-1.mem
printf "OK\n"

printf "[1/2 2/3] Running Scramble BAM-to-CRAM encoder ... "
$time -v -o $input_bam.$crumble_string-1.bam.$scramble_string-1.time $scramble -r $ref_fasta -t $num_threads $input_bam.$crumble_string-1.bam $input_bam.$crumble_string-1.bam.cram &> $input_bam.$crumble_string-1.bam.$scramble_string-1.log
printf "OK\n"

printf "[1/2 3/3] Reporting CRAM size ... "
$cram_size $input_bam.$crumble_string-1.bam.cram &> $input_bam.$crumble_string-1.bam.$scramble_string.cram_size
mv $input_bam.$crumble_string-1.bam.cram $input_bam.$crumble_string-1.bam.$scramble_string
printf "OK\n"

printf "[2/2 1/3] Running Crumble BAM-to-BAM encoder with compression level 9 ... "
$time -v -o $input_bam.$crumble_string-9.time $crumble -v -9 $input_bam $input_bam.$crumble_string-9.bam &> $input_bam.$crumble_string-9.log
printf "OK\n"

printf "[2/2 2/3] Running Scramble BAM-to-CRAM encoder ... "
$time -v -o $input_bam.$crumble_string-9.bam.$scramble_string-1.time $scramble -r $ref_fasta -t $num_threads $input_bam.$crumble_string-9.bam $input_bam.$crumble_string-9.bam.cram &> $input_bam.$crumble_string-9.bam.$scramble_string.log
printf "OK\n"

printf "[2/2 3/3] Reporting CRAM size ... "
$cram_size $input_bam.$crumble_string-9.bam.cram &> $input_bam.$crumble_string-9.bam.$scramble_string.cram_size
mv $input_bam.$crumble_string-9.bam.cram $input_bam.$crumble_string-9.bam.$scramble_string
printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


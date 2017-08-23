#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then printf "Usage: $0 input_bam ref_fasta num_threads\n"; exit -1; fi

input_bam=$1
printf "Input BAM file: $input_bam\n"
ref_fasta=$2
printf "Reference FASTA file: $ref_fasta\n"
num_threads=$3
printf "Number of threads: $num_threads\n"

if [ ! -f $input_bam ]; then printf "Error: Input BAM file $input_bam is not a regular file.\n"; exit -1; fi
if [ ! -f $ref_fasta ]; then printf "Error: Reference FASTA file $ref_fasta is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
cram_size="/project/dna/install/io_lib-1.14.6/bin/cram_size"
crumble="/project/dna/install/crumble-0.5/crumble"
crumble_string="crumble-0.5"
samtools="/project/dna/install/samtools-1.3/bin/samtools"
scramble="/project/dna/install/io_lib-1.14.6/bin/scramble"
scramble_string="scramble-1.14.6"
time="/usr/bin/time"

if [ ! -x $cram_size ]; then printf "Error: Binary file $cram_size is not executable.\n"; exit -1; fi
if [ ! -x $crumble ]; then printf "Error: Binary file $crumble is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
if [ ! -x $scramble ]; then printf "Error: Binary file $scramble is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Constructing FASTA index file\n  from: $ref_fasta\n  to: $ref_fasta.fai\n"
if [ -f $ref_fasta.fai ]; then
    printf "$ref_fasta.fai already exists (not reproducing it)\n"
else
    $samtools faidx $ref_fasta
fi

printf "Crumble BAM-to-BAM encoding with compression level 1\n  from: $input_bam\n  to: $input_bam.$crumble_string-1.bam\n"
cmd="$crumble -v -1 $input_bam $input_bam.$crumble_string-1.bam"
if [ -f $input_bam.$crumble_string-1.bam ]; then
    printf "$input_bam.$crumble_string-1.bam already exists (not reproducing it)\n"
else
    $time -v -o $input_bam.$crumble_string-1.time $cmd &> $input_bam.$crumble_string-1.log
fi

printf "Scramble BAM-to-CRAM encoding\n  from: $input_bam.$crumble_string-1.bam\n  to: $input_bam.$crumble_string-1.bam.cram\n"
if [ -f $input_bam.$crumble_string-1.bam.cram ]; then
    printf "$input_bam.$crumble_string-1.bam.cram already exists (not reproducing it)\n"
else
    $scramble -r $ref_fasta -t $num_threads $input_bam.$crumble_string-1.bam $input_bam.$crumble_string-1.bam.cram &> $input_bam.$crumble_string-1.bam.$scramble_string.log
fi

printf "Reporting CRAM size\n  from: $input_bam.$crumble_string-1.bam.cram\n  to: $input_bam.$crumble_string-1.bam.$scramble_string.cram_size\n"
if [ -f $input_bam.$crumble_string-1.bam.$scramble_string.cram_size ]; then
    printf "$input_bam.$crumble_string-1.bam.$scramble_string.cram_size already exists (not reproducing it)\n"
else
    $cram_size $input_bam.$crumble_string-1.bam.cram &> $input_bam.$crumble_string-1.bam.$scramble_string.cram_size
    mv $input_bam.$crumble_string-1.bam.cram $input_bam.$crumble_string-1.bam.$scramble_string
fi

printf "Crumble BAM-to-BAM encoding with compression level 9\n  from: $input_bam\n  to: $input_bam.$crumble_string-9.bam\n"
cmd="$crumble -v -9 $input_bam $input_bam.$crumble_string-9.bam"
if [ -f $input_bam.$crumble_string-9.bam ]; then
    printf "$input_bam.$crumble_string-9.bam already exists (not reproducing it)\n"
else
    $time -v -o $input_bam.$crumble_string-9.time $cmd &> $input_bam.$crumble_string-9.log
fi

printf "Scramble BAM-to-CRAM encoding\n  from: $input_bam.$crumble_string-9.bam\n  to: $input_bam.$crumble_string-9.bam.cram\n"
if [ -f $input_bam.$crumble_string-9.bam.cram ]; then
    printf "$input_bam.$crumble_string-9.bam.cram already exists (not reproducing it)\n"
else
    $scramble -r $ref_fasta -t $num_threads $input_bam.$crumble_string-9.bam $input_bam.$crumble_string-9.bam.cram &> $input_bam.$crumble_string-9.bam.$scramble_string.log
fi

printf "Reporting CRAM size\n  from: $input_bam.$crumble_string-9.bam.cram\n  to: $input_bam.$crumble_string-9.bam.$scramble_string.cram_size\n"
if [ -f $input_bam.$crumble_string-9.bam.$scramble_string.cram_size ]; then
    printf "$input_bam.$crumble_string-9.bam.$scramble_string.cram_size already exists (not reproducing it)\n"
else
    $cram_size $input_bam.$crumble_string-9.bam.cram &> $input_bam.$crumble_string-9.bam.$scramble_string.cram_size
    mv $input_bam.$crumble_string-9.bam.cram $input_bam.$crumble_string-9.bam.$scramble_string
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


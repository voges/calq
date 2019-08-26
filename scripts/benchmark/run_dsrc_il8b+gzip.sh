#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 input_sam num_threads\n"; exit -1; fi

input_sam=$1
num_threads=$2

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
dsrc="/project/dna/install/dsrc-2.0/dsrc"
dsrc_string="dsrc-2.0_il8b"
gzip="/bin/gzip"
python="/usr/bin/python"
samtools="/project/dna/install/samtools-1.3/bin/samtools"
time="/usr/bin/time"

# Python scripts
sam2fastq_py="/home/muenteferi/Dokumente/calq/src/ngstools/sam2fastq.py"
xtract_part_fastq_py="/home/muenteferi/Dokumente/calq/src/ngstools/xtract_part_fastq.py"
replace_qual_sam_py="/home/muenteferi/Dokumente/calq/src/ngstools/replace_qual_sam.py"



###############################################################################
#                                  Compress                                   #
###############################################################################

printf "SAM-to-FASTQ conversion\n"
$python $sam2fastq_py $input_sam 1>$input_sam.fastq

printf "Compressing with DSRC\n"
$dsrc c -d3 -q2 -b256 -l -t$num_threads $input_sam.fastq $input_sam.fastq.$dsrc_string

printf "Decompressing with DSRC\n"
$dsrc d -t$num_threads $input_sam.fastq.$dsrc_string $input_sam.fastq.$dsrc_string.fastq

printf "Extracting quality values\n"
$python $xtract_part_fastq_py $input_sam.fastq.$dsrc_string.fastq 3 1> $input_sam.fastq.$dsrc_string.fastq.qual

printf "Compressing quality values with gzip\n"
$gzip -9 -c $input_sam.fastq.$dsrc_string.fastq.qual > $input_sam.fastq.$dsrc_string.fastq.qual.gz

printf "Constructing SAM file with reconstructed quality values\n"
$python $replace_qual_sam_py $input_sam $input_sam.fastq.$dsrc_string.fastq.qual 1> $input_sam.$dsrc_string.sam

printf "SAM-to-BAM conversion\n"
$samtools view -bh $input_sam.$dsrc_string.sam > $input_sam.$dsrc_string.bam

printf "BAM index creation\n"
$samtools index -b $input_sam.$dsrc_string.bam $input_sam.$dsrc_string.bai

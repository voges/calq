#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 4 ]; then printf "Usage: $0 input_sam qual_type polyploidy block_size\n"; exit -1; fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
qual_type=$2
printf "Quality value type: $qual_type\n"
polyploidy=$3
printf "Polyploidy: $polyploidy\n"
block_size=$4
printf "Block size: $block_size\n"

if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
calq="/home/muenteferi/Dokumente/calqBuild/calq"
calq_string="calq-haplotype"
python="/usr/bin/python"
samtools="/project/omics/install/samtools-1.3/bin/samtools"
time="/usr/bin/time"

#reference
reference="/data/voges/muenteferi/GATK_bundle-2.8-b37/human_g1k_v37.fasta"

# Python scripts
replace_qual_sam_py="/home/muenteferi/Dokumente/calq/src/ngstools/replace_qual_sam.py"

if [ ! -x $calq ]; then printf "Error: Binary file $calq is not executable.\n"; exit -1; fi
if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi
if [ ! -f $replace_qual_sam_py ]; then printf "Error: Python script $replace_qual_sam_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                          Compress and decompress                            #
###############################################################################

printf "Compressing with CALQ\n"
cmd="$calq -q $qual_type -p $polyploidy -b $block_size $input_sam -o $input_sam.$calq_string -r $reference -f --CALQ-Version v2"
$time -v -o $input_sam.$calq_string.enc.time $cmd &> $input_sam.$calq_string.enc.log

printf "Decompressing with CALQ\n"
cmd="$calq -d -s $input_sam $input_sam.$calq_string -o $input_sam.$calq_string.qual -f"
$time -v -o $input_sam.$calq_string.dec.time $cmd &> $input_sam.$calq_string.dec.log

printf "Constructing SAM file with reconstructed quality values\n"
$python $replace_qual_sam_py $input_sam $input_sam.$calq_string.qual 1> $input_sam.$calq_string.sam

printf "SAM-to-BAM conversion\n"
$samtools view -bh $input_sam.$calq_string.sam > $input_sam.$calq_string.bam

printf "BAM index creation\n"
$samtools index -b $input_sam.$calq_string.bam $input_sam.$calq_string.bai

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


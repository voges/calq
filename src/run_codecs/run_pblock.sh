#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 input_sam p (with p={1,2,4,8,16,32})\n"; exit -1; fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
p=$2
printf "p: $p\n" # p={1,2,4,8,16,32}

if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

prblock_compress="/project/dna/prog/libCSAM-da36a12/CompressQual"
prblock_decompress="/project/dna/prog/libCSAM-da36a12/DecompressQual"
pblock_string="pblock-da36a12"
samtools="/project/dna/install/samtools-1.3/bin/samtools"
time="/usr/bin/time"

if [ ! -x $prblock_compress ]; then printf "Error: Binary file $prblock_compress is not executable.\n"; exit -1; fi
if [ ! -x $prblock_decompress ]; then printf "Error: Binary file $prblock_decompress is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Compressing quality values with p=$p\n"
cmd="$prblock_compress $input_sam -q 1 -m 1 -l $p"
$time -v -o $input_sam.$pblock_string-p$p.enc.time $cmd &> $input_sam.$pblock_string-p$p.enc.log

printf "Decompressing quality values\n"
cmd="$prblock_decompress $input_sam.cqual"
$time -v -o $input_sam.$pblock_string-p$p.dec.time $cmd &> $input_sam.$pblock_string-p$p.dec.log

mv $input_sam.cqual $input_sam.qual.$pblock_string-p$p
mv $input_sam.cqual.qual $input_sam.qual.$pblock_string-p$p.qual

printf "Constructing SAM file with reconstructed quality values\n"
$python $replace_qual_sam_py $input_sam $input_sam.qual.$pblock_string-p$p.qual 1> $input_sam.$pblock_string-p$p.sam

printf "SAM-to-BAM conversion\n"
$samtools view -bh $input_sam.$pblock_string-p$p.sam > $input_sam.$pblock_string-p$p.bam

printf "BAM index creation\n"
$samtools index -b $input_sam.$pblock_string-p$p.bam $input_sam.$pblock_string-p$p.bai

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


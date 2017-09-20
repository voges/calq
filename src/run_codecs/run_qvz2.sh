#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then printf "Usage: $0 input_sam T\n"; exit -1; fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
T=$2
printf "T: $T\n"

if [ ! -f $input_sam ]; then printf "Error: Input SAM file $input_sam is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
python="/usr/bin/python"
qvz2="/project/dna/install/qvz2-d5383c6/qvz2"
qvz2_string="qvz2-d5383c6"
time="/usr/bin/time"

# Python scripts
$replace_qual_sam_py="/home/voges/git/calq/src/ngstools/replace_qual_sam.py"
$xtract_field_sam_py="/home/voges/git/calq/src/ngstools/xtract_field_sam.py"

if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -x $qvz2 ]; then printf "Error: Binary file $qvz2 is not executable.\n"; exit -1; fi
if [ ! -x $time ]; then printf "Error: Binary file $time is not executable.\n"; exit -1; fi
if [ ! -f $replace_qual_sam_py ]; then printf "Error: Python script $replace_qual_sam_py is not a regular file.\n"; exit -1; fi
if [ ! -f $xtract_field_sam_py ]; then printf "Error: Python script $xtract_field_sam_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Extracting quality values from SAM file\n"
$python $xtract_field_sam_py $input_sam 10 1> $input_sam.qual

printf "Running QVZ2\n"
cmd="$qvz2 -t $T -v -u $input_sam.$qvz2_string-t$T.qual $input_sam.qual $input_sam.$qvz2_string-t$T"
$time -v -o $input_sam.$qvz2_string-t$T.time $cmd &> $input_sam.$qvz2_string-t$T.log
wc -c $input_sam.$qvz2_string-t$T &>> $input_sam.$qvz2_string-t$T.log

printf "Constructing new SAM file with QVZ2'd quality values\n"
$python $replace_qual_sam_py $input_sam $input_sam.$qvz2_string-t$T.qual 1>$input_sam.$qvz2_string-t$T.sam

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup ... "
#
printf "OK\n";


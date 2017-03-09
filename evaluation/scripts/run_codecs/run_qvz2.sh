#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 2 ]; then
    printf "Usage: $0 input_sam T\n"
    exit -1
fi

input_sam=$1
printf "Input SAM file: $input_sam\n"
T=$2
printf "T: $T\n"

printf "Checking input SAM file $input_sam ... "
if [ ! -f $input_sam ]; then printf "did not find input SAM file: $input_sam\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

time="/usr/bin/time"
python="/usr/bin/python"
qvz2="/project/dna/install/qvz2-d5383c6/qvz"
qvz2_string="qvz2-d5383c6"
xtract_qual_sam_py="/home/voges/git/ngstools/xtract_qual_sam.py"
replace_qual_sam_py="/home/voges/git/ngstools/replace_qual_sam.py"

printf "Checking executables ... "
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -x $qvz2 ]; then printf "did not find $qvz2\n"; exit -1; fi
if [ ! -e $replace_qual_sam_py ]; then printf "did not find $replace_qual_sam_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

$python $xtract_qual_sam_py $input_sam 2> $input_sam.qual
$time -v -o $input_sam.$qvz2_string-t$t.time $qvz2 -t $T -v -u $input_sam.$qvz2_string-t$T.qual $input_sam.qual $input_sam.$qvz2_string-t$T &> $input_sam.$qvz2_string-t$T.log
wc -c $input_sam.$qvz2_string-t$T &>> $input_sam.$qvz2_string-t$T.log
$python $replace_qual_sam_py $input_sam $input_sam.$qvz2_string-t$T.qual
mv $input_sam.new_qual.sam $input_sam.$qvz2_string-t$T.sam

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


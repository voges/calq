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

# Binaries
pgrep="/usr/bin/pgrep"
python="/usr/bin/python"
qvz2="/project/dna/install/qvz2-d5383c6/qvz2"
qvz2_string="qvz2-d5383c6"
time="/usr/bin/time"

# Python scripts
ps_mem_py="/home/voges/git/calq/evaluation/scripts/ps_mem/ps_mem.py"
replace_qual_sam_py="/home/voges/git/ngstools/replace_qual_sam.py"
xtract_qual_sam_py="/home/voges/git/ngstools/xtract_qual_sam.py"

printf "Checking executables ... "
if [ ! -x $pgrep ]; then printf "did not find $pgrep\n"; exit -1; fi
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -x $qvz2 ]; then printf "did not find $qvz2\n"; exit -1; fi
if [ ! -x $time ]; then printf "did not find $time\n"; exit -1; fi
if [ ! -e $ps_mem_py ]; then printf "did not find $ps_mem_py\n"; exit -1; fi
if [ ! -e $replace_qual_sam_py ]; then printf "did not find $replace_qual_sam_py\n"; exit -1; fi
if [ ! -e $xtract_qual_sam_py ]; then printf "did not find $xtract_qual_sam_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                  Compress                                   #
###############################################################################

printf "Extracting quality values from SAM file ... "
$python $xtract_qual_sam_py $input_sam 2> $input_sam.qual
printf "OK\n"

printf "Running QVZ2 ... "
cmd="$qvz2 -t $T -v -u $input_sam.$qvz2_string-t$T.qual $input_sam.qual $input_sam.$qvz2_string-t$T"
$time -v -o $input_sam.$qvz2_string-t$T.time $cmd &> $input_sam.$qvz2_string-t$T.log &
time_pid=$!
cmd_pid=$($pgrep -P $time_pid)
printf "Command being traced: \"$cmd\"\n" > $input_sam.$qvz2_string-t$T.mem
$python $ps_mem_py -t -w 1 --swap -p $cmd_pid >> $input_sam.$qvz2_string-t$T.mem
wc -c $input_sam.$qvz2_string-t$T &>> $input_sam.$qvz2_string-t$T.log
printf "OK\n"

printf "Constructing new SAM file with QVZ2'd quality values ... "
$python $replace_qual_sam_py $input_sam $input_sam.$qvz2_string-t$T.qual
mv $input_sam.new_qual.sam $input_sam.$qvz2_string-t$T.sam
printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
#
printf "OK\n";


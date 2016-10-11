#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 infile blocksize polyploidy qvtype"
    exit -1
fi

calq="/data/gidb/MPEG/calq"

date;
/usr/bin/time -v -o $1.cq.time $calq -b $2 -f -p $3 -t $4 -v $1 -o $1.cq 2>&1 | tee $1.cq.tee
date;
echo "Wrote timing info to $1.cq.time"
echo "Wrote log to $1.cq.tee"


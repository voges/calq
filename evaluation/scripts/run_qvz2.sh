#!/bin/bash
set -x
/home/voges/git/calq/evaluation/scripts/xtract_qual_sam.py $1 2> $1.qual
/project/dna/prog/qvz2-d5383c6/qvz2 -t $2 -v -u $1.qvz2-d5383c6-t$2.qual $1.qual $1.qvz2-d5383c6-t$2 2>&1 | tee $1.qvz2-d5383c6-t$2.enc.log
/home/voges/git/calq/evaluation/scripts/replace_qual_sam.py $1 $1.qvz2-d5383c6-t$2.qual
mv $1.new_qual.sam $1.qvz2-d5383c6-t$2.sam

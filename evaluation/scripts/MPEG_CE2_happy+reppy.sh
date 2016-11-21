#!/bin/bash

###############################################################################
#                Script to report variant calling results with                #
#                 happy (https://github.com/Illumina/hap.py)                  #
#                                    and                                      #
#              reppy (https://github.com/ga4gh/benchmarking-tools/)           #
###############################################################################

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 num_threads ref.fasta variants.vcf golden.vcf.gz golden.bed"
    exit -1
fi

set -x

num_threads=$1
ref_FASTA=$2
variants_VCF=$3
golden_VCF_GZ=$4
golden_BED=$5

### Programs
install_path="/project/dna/install"
hap_py="$install_path/hap.py-0.3.1/bin/hap.py"
rep_py="$install_path/benchmarking-tools-c458561/reporting/basic/bin/rep.py"

### File names
base=$(echo $variants_VCF | sed 's/\.[^.]*$//') # strip .vcf
happy_root=$base".happy"
reppy_HTML=$base".reppy.html"

### Run hap.py
python $hap_py --threads $num_threads --verbose $golden_VCF_GZ $variants_VCF -f $golden_BED -o $happy_root -r $ref_FASTA --roc VQLSOD

### Run rep.py
rm -f "$happy_root".rep.tsv
printf "method\tcomparisonmethod\tfiles\n" >> "$happy_root".rep.tsv
printf "$variants_VCF\t$golden_VCF_GZ\t" >> "$happy_root".rep.tsv
printf "$happy_root.extended.csv," >> "$happy_root".rep.tsv
for i in "$happy_root".roc.Locations.*; do printf "$i,"; done >> "$happy_root".rep.tsv
sed -i '$ s/.$//' "$happy_root".rep.tsv
python $rep_py -o $reppy_HTML -l "$happy_root".rep.tsv


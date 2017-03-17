#!/bin/bash

###############################################################################
#                Script to report variant calling results with                #
#                hap.py (https://github.com/Illumina/hap.py)                  #
#                                    and                                      #
#             rep.py (https://github.com/ga4gh/benchmarking-tools/)           #
###############################################################################

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 6 ]; then
    printf "Usage: $0 num_threads log_txt variants_vcf ref_fasta golden_vcf_gz golden_bed\n"
    exit -1
fi

num_threads=$1
printf "Number of threads: $num_threads\n"
log_txt=$2
printf "Log file: $log_txt\n"
variants_vcf=$3
printf "Input VCF file: $variants_vcf\n"
ref_fasta=$4
printf "Reference FASTA file: $ref_fasta\n"
golden_vcf_gz=$5
printf "Golden VCF file (gzipped): $golden_vcf_gz\n"
golden_bed=$6
printf "Golden BED file: $golden_bed\n"

printf "Checking log file $log_txt ... "
if [ -f $log_txt ]; then printf "log file already exists: $log_txt\n"; exit -1; fi
printf "OK\n"

printf "Checking input VCF file $variants_vcf ... "
if [ ! -f $variants_vcf ]; then printf "did not find input VCF file: $variants_vcf\n"; exit -1; fi
printf "OK\n"

printf "Checking reference FASTA file $ref_fasta ... "
if [ ! -f $ref_fasta ]; then printf "did not find reference FASTA file: $ref_fasta\n"; exit -1; fi
printf "OK\n"

printf "Checking golden VCF file (gzipped) $golden_vcf_gz ... "
if [ ! -f $golden_vcf_gz ]; then printf "did not find golden VCF file (gzipped): $golden_vcf_gz\n"; exit -1; fi
printf "OK\n"

printf "Checking golden BED file $golden_bed ... "
if [ ! -f $golden_bed ]; then printf "did not find golden BED file: $golden_bed\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

python="/usr/bin/python"
hap_py="/project/dna/install/hap.py-0.3.1/bin/hap.py"
rep_py="/project/dna/install/benchmarking-tools-c458561/reporting/basic/bin/rep.py"

printf "Checking executables ... "
if [ ! -x $python ]; then printf "did not find $python\n"; exit -1; fi
if [ ! -e $hap_py ]; then printf "did not find $hap_py\n"; exit -1; fi
if [ ! -e $rep_py ]; then printf "did not find $rep_py\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                 Report results with hap.py and rep.py                       #
###############################################################################

printf "[1/2] Running hap.py ... "
$python $hap_py --threads $num_threads --verbose $golden_vcf_gz $variants_vcf -f $golden_bed -o $variants_vcf.happy -r $ref_fasta --roc VQLSOD &>>$log_txt
printf "OK\n"

printf "[2/2] Running rep.py ... "
rm -f "$variants_vcf".happy.rep.tsv
printf "method\tcomparisonmethod\tfiles\n" >> "$variants_vcf".happy.rep.tsv
printf "$variants_vcf\t$golden_vcf_gz\t" >> "$variants_vcf".happy.rep.tsv
printf "$variants_vcf.happy.extended.csv," >> "$variants_vcf".happy.rep.tsv
for i in "$variants_vcf".happy.roc.Locations.*; do printf "$i,"; done >> "$variants_vcf".happy.rep.tsv
sed -i '$ s/.$//' "$variants_vcf".happy.rep.tsv
$python $rep_py -o "$variants_vcf".happy.rep.tsv.reppy.html -l "$variants_vcf".happy.rep.tsv &>>$log_txt
printf "OK\n"

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
rm -f $variants_vcf.happy.extended.csv
rm -f $variants_vcf.happy.metrics.json
rm -f $variants_vcf.happy.roc.*
rm -f $variants_vcf.happy.vcf*
printf "OK\n";


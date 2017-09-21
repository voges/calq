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

if [ -f $log_txt ]; then printf "Error: File $log_txt file already exists.\n"; exit -1; fi
if [ ! -f $variants_vcf ]; then printf "Error: Input VCF file $variants_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $ref_fasta ]; then print "Error: Reference FASTA file $ref_fasta is not a regular file.\n"; exit -1; fi
if [ ! -f $golden_vcf_gz ]; then printf "Error: Golden VCF file (gzipped) $golden_vcf_gz is not a regular file.\n"; exit -1; fi
if [ ! -f $golden_bed ]; then printf "Error: Golden BED file $golden_bed is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
python="/usr/bin/python"

# Python scripts
hap_py="/project/dna/install/hap.py-0.3.1/bin/hap.py"
rep_py="/project/dna/install/benchmarking-tools-c458561/reporting/basic/bin/rep.py"

if [ ! -x $python ]; then printf "Error: Binary file $python is not executable.\n"; exit -1; fi
if [ ! -f $hap_py ]; then printf "Error: Python script $hap_py is not a regular file.\n"; exit -1; fi
if [ ! -f $rep_py ]; then printf "Error: Python script $rep_py is not a regular file.\n"; exit -1; fi

###############################################################################
#                 Report results with hap.py and rep.py                       #
###############################################################################

printf "[1/2] Running hap.py\n"
$python $hap_py --threads $num_threads --verbose $golden_vcf_gz $variants_vcf -f $golden_bed -o $variants_vcf.happy -r $ref_fasta --roc VQLSOD &>>$log_txt

printf "[2/2] Running rep.py\n"
rm -f "$variants_vcf".happy.rep.tsv
printf "method\tcomparisonmethod\tfiles\n" >> "$variants_vcf".happy.rep.tsv
printf "$variants_vcf\t$golden_vcf_gz\t" >> "$variants_vcf".happy.rep.tsv
printf "$variants_vcf.happy.extended.csv," >> "$variants_vcf".happy.rep.tsv
for i in "$variants_vcf".happy.roc.Locations.*; do printf "$i,"; done >> "$variants_vcf".happy.rep.tsv
sed -i '$ s/.$//' "$variants_vcf".happy.rep.tsv
$python $rep_py -o "$variants_vcf".happy.rep.tsv.reppy.html -l "$variants_vcf".happy.rep.tsv &>>$log_txt

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup\n"
rm -f $variants_vcf.happy.extended.csv
rm -f $variants_vcf.happy.metrics.json
rm -f $variants_vcf.happy.roc.*
rm -f $variants_vcf.happy.vcf*
printf "Done\n";


#!/bin/bash

###############################################################################
#                     Extract chromosome from VCF file                        #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 file.vcf ref.fa chromosome"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .vcf
ref=$2
chromosome=$3

install_path="/project/dna/install"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"

date; echo "Extracting chromosome $chromosome from VCF file $1"

java -jar $GenomeAnalysisTK_jar \
   -T SelectVariants \
   -R $ref \
   -V $root.vcf \
   -o $root.$chromosome.vcf \
   -L $chromosome

date; echo "Wrote variants on chromosome $chromosome to VCF file $root.$chromosome.vcf"


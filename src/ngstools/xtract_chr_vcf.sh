#!/bin/bash

###############################################################################
#                     Extract chromosome from VCF file                        #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 file.vcf ref.[fa|fasta] chromosome"
    exit -1
fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .vcf
ref=$2
chromosome=$3

install_path="/project/dna/install"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"
samtools="$install_path/samtools-1.3/bin/samtools"

if [ -f $2.fai ]; then
    date; echo "FASTA index file $2.fai already exists. Not reproducing it."
    #touch $2.fai
else
    date; echo "Constructing FASTA index file $2.fai"
    $samtools faidx $2
fi

date; echo "Extracting chromosome $chromosome from VCF file $1"

java -jar $GenomeAnalysisTK_jar \
   -T SelectVariants \
   -R $ref \
   -V $root.vcf \
   -o $root.$chromosome.vcf \
   -L $chromosome

if [ -f $root.$chromosome.vcf ]; then
    date; echo "Wrote variants on chromosome $chromosome to VCF file $root.$chromosome.vcf"
else
    date; echo "Extracting failed"
fi

#!/bin/bash

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

echo "Creating FASTA index file $2.fai"
if [ -f $2.fai ]; then
    echo "FASTA index file $2.fai already exists. Not reproducing it."
else
    $samtools faidx $2
fi

echo "Writing chromosome $chromosome from VCF file $1 to VCF file $root.$chromosome.vcf"
if [ -f $root.$chromosome.vcf ]; then
    echo "VCF file $root.$chromosome.vcf already exists. Not reproducing it."
else
    java -jar $GenomeAnalysisTK_jar \
       -T SelectVariants \
       -R $ref \
       -V $root.vcf \
       -o $root.$chromosome.vcf \
       -L $chromosome
        if [ -f $root.$chromosome.vcf ]; then
    else
        echo "Writing failed"
    fi
fi


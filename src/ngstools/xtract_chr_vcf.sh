#!/bin/bash

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then printf "Usage: $0 file.vcf ref.[fa|fasta] chromosome\n"; exit -1; fi

root=$(echo $1 | sed 's/\.[^.]*$//') # strip .vcf
ref=$2
chromosome=$3

printf "Input VCF file: $1\n"
printf "Reference FASTA file: $ref\n"
printf "Chromosome: $chromosome\n"

if [ ! -f $1 ]; then printf "Error: Input VCF file $1 is not a regular file.\n"; exit -1; fi
if [ ! -f $ref ]; then printf "Error: Reference FASTA file $ref is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

install_path="/project/dna/install"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"
java="/usr/bin/java"
java_opts=""
samtools="$install_path/samtools-1.3/bin/samtools"

if [ ! -x $GenomeAnalysisTK_jar ]; then printf "Error: JAR file $GenomeAnalysisTK_jar is not executable.\n"; exit -1; fi
if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi

###############################################################################
#                                 Extraction                                  #
###############################################################################

printf "Creating FASTA index file: $ref.fai\n"
if [ -f $ref.fai ]; then
    printf "$ref.fai already exists (not reproducing it)\n"
else
    $samtools faidx $2
fi

printf "Extracting\n  chromosome: $chromosome\n  from: $1\n  to: $root.$chromosome.vcf\n"
if [ -f $root.$chromosome.vcf ]; then
    printf "$root.$chromosome.vcf already exists (not reproducing it)\n"
else
    $java $java_opts -jar $GenomeAnalysisTK_jar \
        -T SelectVariants \
        -R $ref \
        -V $root.vcf \
        -o $root.$chromosome.vcf \
        -L $chromosome
fi

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#
printf "Done\n"


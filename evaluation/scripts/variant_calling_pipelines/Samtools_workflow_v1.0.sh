#!/bin/bash

###############################################################################
#    Script for performing WGS/WES Mapping to Variant Calls - Version 1.0     #
#                      http://www.htslib.org/workflow/                        #
###############################################################################

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 num_threads input_fastq sample"
    exit -1
fi

num_threads=$1
input_fastq=$2
sample=$3
root=$(echo $input_fastq | sed 's/\.[^.]*$//') # strip .fastq

###############################################################################
#                          Data and programs                                  #
###############################################################################

# GATK bundle
gatk_bundle_path="/data/voges/MPEG/GATK_bundle-2.8-b37"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
mills_vcf="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"

# Binaries
bwa="/project/dna/install/bwa-0.7.13/bwa"
java="/usr/bin/java"
samtools="/project/dna/install/samtools-1.3/bin/samtools"
tabix="/project/dna/install/project/dna/install/htslib-1.3/bin/tabix"

# JAR files
picard_jar="/project/dna/install/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"

# Reference indexing
if [ ! -e "${ref_fasta}.fai" ]; then
    $samtools faidx $ref_fasta
fi

###############################################################################
#                                  Mapping                                    #
###############################################################################

# Generate BWA index
$bwa index -a bwtsw $ref_fasta

# BWA MEM alignment
$bwa mem -t $num_threads -M $ref_fasta $input_fastq > $root.aln_bwa.sam

# Clean up read pairing information and flags and convert to BAM
$samtools fixmate -O bam $root.aln_bwa.sam $root.aln_bwa.fixmate.bam
rm -f $root.aln_bwa.sam

# Sort BAM file
$samtools sort -@ $num_threads -O bam $root.aln_bwa.fixmate.bam > $root.aln_bwa.fixmate.sorted.bam
rm -f $root.aln_bwa.fixmate.bam

# Add read group name
$java -jar $picard_jar AddOrReplaceReadGroups I=$root.aln_bwa.fixmate.sorted.bam O=$root.aln_bwa.fixmate.sorted.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample
rm -f $root.aln_bwa.fixmate.sorted.bam

###############################################################################
#                                 Improvement                                 #
###############################################################################

# Reduce the number of miscalls of INDELs by realigning
$java -Xmx2g -jar $GenomeAnalysisTK_jar -T RealignerTargetCreator -R $ref_fasta -I $root.aln_bwa.fixmate.sorted.rg.bam -o $root.realigner_target.intervals --known $mills_vcf
$java -Xmx4g -jar $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_fasta -I $root.aln_bwa.fixmate.sorted.rg.bam -targetIntervals $root.realigner_target.intervals -o $root.aln_bwa.fixmate.sorted.rg.realn.bam
rm -f $root.realigner_target.intervals
rm -f $root.aln_bwa.fixmate.sorted.rg.bam

# Mark duplicates in the BAM file
$java -Xmx2g -jar $picard_jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=$root.aln_bwa.fixmate.sorted.rg.realn.bam O=$root.aln_bwa.fixmate.sorted.rg.realn.dupmark.bam M=$root.dedup_metrics.txt
rm -f $root.dedup_metrics.txt
rm -f $root.aln_bwa.fixmate.sorted.rg.realn.bam

### Index the BAM file
$samtools index $root.aln_bwa.fixmate.sorted.rg.realn.dupmark.bam

###############################################################################
#                               Variant calling                               #
###############################################################################

# Call variants
$samtools mpileup -ugf $ref_fasta $root.aln_bwa.fixmate.sorted.rg.realn.dupmark.bam | $bcftools call -vmO z -o $root.raw_variants.vcf.gz

# Index VCF file
$tabix -p vcf $root.raw_variants.vcf.gz


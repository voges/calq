#!/bin/bash

###############################################################################
#    Script for performing WGS/WES Mapping to Variant Calls - Version 1.0     #
#                      http://www.htslib.org/workflow/                        #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 num_threads reads.fastq sample"
    exit -1
fi

###############################################################################
#                               Command line                                  #
###############################################################################
date
set -x
script_name=$0
num_threads=$1
reads_FASTQ=$2
sample=$3
root=$(echo $reads_FASTQ | sed 's/\.[^.]*$//') # strip .fastq
root="$root.$script_name-$sample"

###############################################################################
#                          Data and programs                                  #
###############################################################################

### GATK bundle
gatk_bundle_path="/data/gidb/GATK_bundle-2.8-b37"
ref_FASTA="$gatk_bundle_path/human_g1k_v37.fasta"
mills_VCF="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"

### Programs
install_path="/project/dna/install"
bwa="$install_path/bwa-0.7.13/bwa"
samtools="$install_path/samtools-1.3/bin/samtools"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"
tabix="$install_path/project/dna/install/htslib-1.3/bin/tabix"

### Temporary directories
javaIOTmpDir="$root/javaIOTmp.dir/"
samtoolsTmpDir="$root/samtoolsTmp.dir/"

### Reference indexing
if [ ! -e "${ref_FASTA}.fai" ]; then
    $samtools faidx $ref_FASTA; date
fi

###############################################################################
#                                  Mapping                                    #
###############################################################################

### Generate BWA index
date; $bwa index -a bwtsw $ref_FASTA

### BWA MEM alignment
date; $bwa mem -t $num_threads -M $ref_FASTA $reads_FASTQ > $root.aln_bwa.sam

### Clean up read pairing information and flags and convert to BAM
date; $samtools fixmate -O bam $root.aln_bwa.sam $root.aln_bwa.fixmate.bam
rm -f $root.aln_bwa.sam

### Sort BAM file
date; $samtools sort -T $samtoolsTmpDir -@ $num_threads -O bam $root.aln_bwa.fixmate.bam > $root.aln_bwa.fixmate.sorted.bam
rm -f $root.aln_bwa.fixmate.bam

### Add read group name
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$root.aln_bwa.fixmate.sorted.bam O=$root.aln_bwa.fixmate.sorted.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample
rm -f $root.aln_bwa.fixmate.sorted.bam

###############################################################################
#                                 Improvement                                 #
###############################################################################

### Reduce the number of miscalls of INDELs by realigning
date; java -Xmx2g -jar $GenomeAnalysisTK_jar -T RealignerTargetCreator -R $ref_FASTA -I $root.aln_bwa.fixmate.sorted.rg.bam -o $root.realigner_target.intervals --known $mills_VCF
date; java -Xmx4g -jar $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_FASTA -I $root.aln_bwa.fixmate.sorted.rg.bam -targetIntervals $root.realigner_target.intervals -o $root.aln_bwa.fixmate.sorted.rg.realn.bam
rm -f $root.realigner_target.intervals
rm -f $root.aln_bwa.fixmate.sorted.rg.bam

### Mark duplicates in the BAM file
date; java -Xmx2g -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=$root.aln_bwa.fixmate.sorted.rg.realn.bam O=$root.aln_bwa.fixmate.sorted.rg.realn.dupmark.bam M=$root.dedup_metrics.txt
rm -f $root.dedup_metrics.txt
rm -f $root.aln_bwa.fixmate.sorted.rg.realn.bam

### Index the BAM file
date; $samtools index $root.aln_bwa.fixmate.sorted.rg.realn.dupmark.bam

###############################################################################
#                               Variant calling                               #
###############################################################################

### Call variants
date; $samtools mpileup -ugf $ref_FASTA $root.aln_bwa.fixmate.sorted.rg.realn.dupmark.bam | $bcftools call -vmO z -o $root.raw_variants.vcf.gz

### Index VCF file
date; $tabix -p vcf $root.raw_variants.vcf.gz

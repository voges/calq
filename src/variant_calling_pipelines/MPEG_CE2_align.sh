#!/bin/bash

###############################################################################
#   Script for performing the ISO/IEC JTC 1/SC 29/WG 11 and ISO/TC 276/WG 5   #
#   Benchmark framework for lossy compression of sequencing quality values    #
#                       (document no. N16727/N156)                            #
#                                    -                                        #
#                                Alignment                                    #
###############################################################################

if [ "$#" -ne 4 ]; then printf "Usage: $0 num_threads reads pairing platform\n"; exit -1; fi

###############################################################################
#                               Command line                                  #
###############################################################################

num_threads=$1
reads=$2
pairing=$3 # 'paired' or 'unpaired'
platform=$4 # ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
root="$reads"

###############################################################################
#                          Data and programs                                  #
###############################################################################

# GATK bundle
gatk_bundle_path="/phys/intern2/MPEG/GATK_bundle-2.8-b37"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_vcf="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_vcf="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
KG_vcf="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_vcf="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_vcf="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

# Binaries
bowtie2="/project/dna/install/bowtie2-2.2.5/bowtie2"
java="/usr/bin/java"
samtools="/project/dna/install/samtools-1.3/bin/samtools"

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
picard_jar="/project/dna/install/picard-tools-2.4.1/picard.jar"

###############################################################################
#                           Alignment with Bowtie2                            #
###############################################################################

$bowtie2-build $ref_fasta $root.bowtie2_idx
if [ "$pairing" = "unpaired" ]; then
    $bowtie2 -x $root.bowtie2_idx -U $reads.fastq -S $root.aln_bowtie2.sam --threads $num_threads
else
    if [ "$pairing" = "paired" ]; then
        $bowtie2 -x $root.bowtie2_idx -1 $reads\_1.fastq -2 $reads\_2.fastq -S $root.aln_bowtie2.sam --threads $num_threads
    else
        echo "pairing argument must be either 'unpaired' or 'paired'"
        exit -1
    fi
fi
#rm -f $root.bowtie2_idx*

###############################################################################
#                             Sorting & indexing                              #
###############################################################################

# Convert SAM to BAM
$samtools view -@ $num_threads -bh $root.aln_bowtie2.sam > $root.aln_bowtie2.bam

# Sort and index BAM file
$samtools sort -@ $num_threads -O bam $root.aln_bowtie2.bam > $root.aln_bowtie2.sorted.bam
$samtools index $root.aln_bowtie2.sorted.bam

###############################################################################
#                              Duplicate removal                              #
###############################################################################

# Mark duplicates in the BAM file
$java -jar $picard_jar MarkDuplicates I=$root.aln_bowtie2.sorted.bam O=$root.aln_bowtie2.sorted.dupmark.bam M=$root.dedup_metrics.txt ASSUME_SORTED=true

# Label the BAM headers and index the resulting file
$java -jar $picard_jar AddOrReplaceReadGroups I=$root.aln_bowtie2.sorted.dupmark.bam O=$root.aln_bowtie2.sorted.dupmark.rg.bam RGID=1 RGLB=Library RGPL=$platform RGPU=PlatformUnit RGSM=$root
$samtools index $root.aln_bowtie2.sorted.dupmark.rg.bam

###############################################################################
#                              Indel realignment                              #
###############################################################################

$java -jar $GenomeAnalysisTK_jar -T RealignerTargetCreator -nt $num_threads -R $ref_fasta -I $root.aln_bowtie2.sorted.dupmark.rg.bam --known $mills_vcf -o $root.realigner_target.intervals
$java -jar $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_fasta -I $root.aln_bowtie2.sorted.dupmark.rg.bam -targetIntervals $root.realigner_target.intervals -o $root.aln_bowtie2.sorted.dupmark.rg.realn.bam

###############################################################################
#                      Base quality score recalibration                       #
###############################################################################

$java -jar $GenomeAnalysisTK_jar -T BaseRecalibrator -nct $num_threads -R $ref_fasta -I $root.aln_bowtie2.sorted.dupmark.rg.realn.bam -knownSites $dbsnps_vcf -knownSites $mills_vcf -knownSites $indels_vcf -o $root.bqsr.table
$java -jar $GenomeAnalysisTK_jar -T PrintReads -nct $num_threads -R $ref_fasta -I $root.aln_bowtie2.sorted.dupmark.rg.realn.bam -BQSR $root.bqsr.table -o $root.aln_bowtie2.sorted.dupmark.rg.realn.recal.bam


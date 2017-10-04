#!/bin/bash

###############################################################################
#   Script for performing the ISO/IEC JTC 1/SC 29/WG 11 and ISO/TC 276/WG 5   #
#   Benchmark framework for lossy compression of sequencing quality values    #
#                       (document no. N16727/N156)                            #
#                                    -                                        #
#                                Alignment                                    #
###############################################################################

if [ "$#" -ne 4 ]; then
    printf "Usage: $0 num_threads reads_root pairing platform\n";
    exit -1;
fi

###############################################################################
#                               Command line                                  #
###############################################################################

num_threads=$1
reads_root=$2
pairing=$3 # 'paired' or 'unpaired'
platform=$4 # ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/phys/ssd/voges/MPEG/GATK_bundle-2.8-b37"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_vcf="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_vcf="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
KG_vcf="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_vcf="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_vcf="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

if [ ! -f $ref_fasta ]; then printf "Error: File $ref_fasta is not a regular file.\n"; exit -1; fi
if [ ! -f $hapmap_vcf ]; then printf "Error: File $hapmap_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $omni_vcf ]; then printf "Error: File $omni_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $KG_vcf ]; then printf "Error: File $KG_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $dbsnps_vcf ]; then printf "Error: File $dbsnps_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $mills_vcf ]; then printf "Error: File $mills_vcf is not a regular file.\n"; exit -1; fi
if [ ! -f $indels_vcf ]; then printf "Error: File $indels_vcf is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
bowtie2="/project/dna/install/bowtie2-2.2.5/bowtie2"
java="/usr/bin/java"
samtools="/project/dna/install/samtools-1.3/bin/samtools"

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
picard_jar="/project/dna/install/picard-tools-2.4.1/picard.jar"

if [ ! -x $bowtie2 ]; then printf "Error: Binary file $bowtie2 is not executable.\n"; exit -1; fi
if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
if [ ! -f $GenomeAnalysisTK_jar ]; then printf "Error: JAR file $GenomeAnalysisTK_jar is not a regular file.\n"; exit -1; fi
if [ ! -f $picard_jar ]; then printf "Error: JAR file $picard_jar is not a regular file.\n"; exit -1; fi

###############################################################################
#                           Alignment with Bowtie2                            #
###############################################################################

$bowtie2-build $ref_fasta $reads_root.bowtie2_idx
if [ "$pairing" = "unpaired" ]; then
    $bowtie2 -x $reads_root.bowtie2_idx -U $reads_root.fastq -S $reads_root.aln_bowtie2.sam --threads $num_threads
else
    if [ "$pairing" = "paired" ]; then
        $bowtie2 -x $reads_root.bowtie2_idx -1 $reads_root\_1.fastq -2 $reads_root\_2.fastq -S $reads_root.aln_bowtie2.sam --threads $num_threads
    else
        echo "pairing argument must be either 'unpaired' or 'paired'"
        exit -1
    fi
fi

###############################################################################
#                             Sorting & indexing                              #
###############################################################################

# Convert SAM to BAM and index BAM file
$samtools view -@ $num_threads -bh $reads_root.aln_bowtie2.sam > $reads_root.aln_bowtie2.bam
$samtools index -b $reads_root.aln_bowtie2.sorted.bam $reads_root.aln_bowtie2.sorted.bai

# Sort and index BAM file
$samtools sort -@ $num_threads -O bam $reads_root.aln_bowtie2.bam > $reads_root.aln_bowtie2.sorted.bam
$samtools index -b $reads_root.aln_bowtie2.sorted.bam $reads_root.aln_bowtie2.sorted.bai

###############################################################################
#                              Duplicate removal                              #
###############################################################################

# Mark duplicates in the BAM file and index the new BAM file
$java -jar $picard_jar MarkDuplicates I=$reads_root.aln_bowtie2.sorted.bam O=$reads_root.aln_bowtie2.sorted.dupmark.bam M=$reads_root.dedup_metrics.txt ASSUME_SORTED=true
$samtools index -b $reads_root.aln_bowtie2.sorted.dupmark.bam $reads_root.aln_bowtie2.sorted.dupmark.bai

# Label the BAM headers and index the resulting file
$java -jar $picard_jar AddOrReplaceReadGroups I=$reads_root.aln_bowtie2.sorted.dupmark.bam O=$reads_root.aln_bowtie2.sorted.dupmark.rg.bam RGID=1 RGLB=Library RGPL=$platform RGPU=PlatformUnit RGSM=$reads_root
$samtools index -b $reads_root.aln_bowtie2.sorted.dupmark.rg.bam $reads_root.aln_bowtie2.sorted.dupmark.rg.bai

###############################################################################
#                              Indel realignment                              #
###############################################################################

$java -jar $GenomeAnalysisTK_jar -T RealignerTargetCreator -nt $num_threads -R $ref_fasta -I $reads_root.aln_bowtie2.sorted.dupmark.rg.bam --known $mills_vcf -o $reads_root.realigner_target.intervals
$java -jar $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_fasta -I $reads_root.aln_bowtie2.sorted.dupmark.rg.bam -targetIntervals $reads_root.realigner_target.intervals -o $reads_root.aln_bowtie2.sorted.dupmark.rg.realn.bam

###############################################################################
#                      Base quality score recalibration                       #
###############################################################################

$java -jar $GenomeAnalysisTK_jar -T BaseRecalibrator -nct $num_threads -R $ref_fasta -I $reads_root.aln_bowtie2.sorted.dupmark.rg.realn.bam -knownSites $dbsnps_vcf -knownSites $mills_vcf -knownSites $indels_vcf -o $reads_root.bqsr.table
$java -jar $GenomeAnalysisTK_jar -T PrintReads -nct $num_threads -R $ref_fasta -I $reads_root.aln_bowtie2.sorted.dupmark.rg.realn.bam -BQSR $reads_root.bqsr.table -o $reads_root.aln_bowtie2.sorted.dupmark.rg.realn.recal.bam


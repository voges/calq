#!/bin/bash

###############################################################################
#   Script for performing the GATK Best Practices for Germline SNP & Indel    #
#               Discovery in Whole Genome and Exome Sequence                  #
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

if [ ! -f $input_fastq ]; then printf "Error: Input FASTQ file $input_fastq is not a regular file.\n"; exit -1; fi

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/phys/intern2/muenteferi/GATK"
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
bwa="/project/dna/install/bwa-0.7.13/bwa"
java="/usr/bin/java"
samtools="/project/dna/install/samtools-1.3/bin/samtools"

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
picard_jar="/project/dna/install/picard-tools-2.4.1/picard.jar"

if [ ! -x $bwa ]; then printf "Error: Binary file $bwa is not executable.\n"; exit -1; fi
if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "Error: Binary file $samtools is not executable.\n"; exit -1; fi
if [ ! -f $GenomeAnalysisTK_jar ]; then printf "Error: JAR file $GenomeAnalysisTK_jar is not a regular file.\n"; exit -1; fi
if [ ! -f $picard_jar ]; then printf "Error: JAR file $picard_jar is not a regular file.\n"; exit -1; fi

###############################################################################
#                       Map and mark duplicates                               #
###############################################################################

# Reference indexing
if [ ! -f "${ref_fasta}.fai" ]; then
    $samtools faidx $ref_fasta
fi

# Generate BWA index
$bwa index -a bwtsw $ref_fasta

# BWA MEM alignment
$bwa mem -t $num_threads -M $ref_fasta $input_fastq > $root.aln_bwa.sam

# Sort SAM file and convert to BAM
$java -jar $picard_jar SortSam I=$root.aln_bwa.sam O=$root.aln_bwa.sorted.bam SORT_ORDER=coordinate

# Mark duplicates in the BAM file
$java -jar $picard_jar MarkDuplicates I=$root.aln_bwa.sorted.bam O=$root.aln_bwa.sorted.dupmark.bam M=$root.dedup_metrics.txt ASSUME_SORTED=true

# Add read group name
$java -jar $picard_jar AddOrReplaceReadGroups I=$root.aln_bwa.sorted.dupmark.bam O=$root.aln_bwa.sorted.dupmark.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample

# Index the BAM file
$java -jar $picard_jar BuildBamIndex I=$root.aln_bwa.sorted.dupmark.rg.bam

###############################################################################
#                      Recalibrate base quality scores                        #
###############################################################################

#$java -jar $GenomeAnalysisTK_jar -nct $num_threads -T baseRecalibrator -R $ref_fasta -I $root.aln_bwa.sorted.dupmark.rg.bam -knownSites $dbsnps_vcf -knownSites $mills_vcf -knownSites $indels_vcf -o $root.bqsr.table
#$java -jar $GenomeAnalysisTK_jar -nct $num_threads -T PrintReads -R $ref_fasta -I $root.aln_bwa.sorted.dupmark.rg.bam -BQSR $root.bqsr.table -o $root.aln_bwa.sorted.dupmark.rg.recal.bam

###############################################################################
#                      Call variants using Haplotype Caller                   #
###############################################################################

SEC=10
SCC=30
$java -jar $GenomeAnalysisTK_jar -T HaplotypeCaller -R $ref_fasta -I $root.aln_bwa.sorted.dupmark.rg.bam --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $root.raw_variants.vcf

###############################################################################
#                      Recalibrate variant quality scores                     #
###############################################################################

resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_vcf"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $KG_vcf"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_vcf"
# Added -mG 4 (default: 8) and -minNumBad 5000 (default: 1000)
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mG 4 -minNumBad 5000"
resourceIndels="mills,known=true,training=true,truth=true,prior=12.0 $mills_vcf"
recalParamsIndels="-an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum"
filterLevel="99.0"

# Argument -tranche is MPEG's theta
$java -jar $GenomeAnalysisTK_jar -R $ref_fasta -T VariantRecalibrator -input $root.raw_variants.vcf -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $root.snps.recal -tranchesFile $root.snps.tranches -rscriptFile $root.snps.r
$java -jar $GenomeAnalysisTK_jar -R $ref_fasta -T ApplyRecalibration -input $root.raw_variants.vcf -mode SNP -recalFile $root.snps.recal -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel -o $root.recalibrated_snps+raw_indels.vcf

# Argument -tranche is MPEG's theta
$java -jar $GenomeAnalysisTK_jar -R $ref_fasta -T VariantRecalibrator -input $root.recalibrated_snps+raw_indels.vcf -resource:$resourceIndels $recalParamsIndels -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $root.indels.recal -tranchesFile $root.indels.tranches -rscriptFile $root.indels.r
$java -jar $GenomeAnalysisTK_jar -R $ref_fasta -T ApplyRecalibration -input $root.recalibrated_snps+raw_indels.vcf -mode INDEL -recalFile $root.indels.recal -tranchesFile $root.indels.tranches --ts_filter_level $filterLevel -o $root.recalibrated_variants.vcf


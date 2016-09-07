#!/bin/bash

###############################################################################
#     Script for performing GATK Best Practices for Germline SNP & Indel      #
#               Discovery in Whole Genome and Exome Sequence                  #
#     https://software.broadinstitute.org/gatk/best-practices/index.php       #
###############################################################################

### More info: https://software.broadinstitute.org/gatk/events/slides/1506/GATKwr8-A-3-GATK_Best_Practices_and_Broad_pipelines.pdf

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 num_threads reads.fastq"
    exit -1
fi

###############################################################################
#                               Command line                                  #
###############################################################################
set +x;echo "";echo "### Command line ###";echo "";set -x

num_threads=$1
reads_FASTQ=$2

base=$(echo $reads_FASTQ | sed 's/\.[^.]*$//') # strip .fastq

###############################################################################
#                                  Data                                       #
###############################################################################
set +x;echo "";echo "### Data ###";echo "";set -x

### GATK bundle
gatk_bundle_path="/data/gidb/simulations/GATK_bundle-2.8-b37"
ref_FASTA="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_VCF="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_VCF="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
snps_VCF="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_VCF="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_VCF="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_VCF="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

### File names
aln_SAM="$base.aln.sam"
aln_sorted_BAM="$base.aln.sorted.bam"
aln_sorted_dedup_BAM="$base.aln.sorted.dedup.bam"
aln_sorted_dedup_rg_BAM="$base.aln.sorted.dedup.rg.bam"
dedupMetrics_TXT="$base.dedupMetrics.txt"
#aln_sorted_dedup_rg_recal_BAM="$base.aln.sorted.dedup.rg.recal.bam"
#bqsr_TABLE="$base.bqsr.table"
raw_variants_VCF="$base.raw_variants.vcf"
snps_RECAL="$base.snps.recal"
snps_TRANCHES="$base.snps.tranches"
snps_R="$base.snps.R"
indels_RECAL="$base.indels.recal"
indels_TRANCHES="$base.indels.tranches"
indels_R="$base.indels.R"
recalibrated_snps_raw_indels_VCF="$base.recalibrated_snps_raw_indels.vcf"
recalibrated_variants_VCF="$base.recalibrated_variants.vcf"

###############################################################################
#                                Programs                                     #
###############################################################################
set +x;echo "";echo "### Programs ###";echo "";set -x

install_path="/project/dna/install"
bwa="$install_path/bwa-0.7.13/bwa"
samtools="$install_path/samtools-1.3/bin/samtools"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"

javaIOTmpDir="$base/javaIOTmp.dir"
samtoolsTmpDir="$basesamtoolTmp.dir"

###############################################################################
#                       Map and mark duplicates                               #
###############################################################################
set +x;echo "";echo "### Map and mark duplicates ###";echo "";set -x

### Generate BWA index
$bwa index -a bwtsw $ref_FASTA

### BWA MEM alignment
$bwa mem -t $num_threads -M $ref_FASTA $reads_FASTQ > $aln_SAM

### Sort SAM file and convert to BAM
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar SortSam I=$aln_SAM O=$aln_sorted_BAM SORT_ORDER=coordinate
rm -f $aln_SAM

### Mark duplicates in the BAM file
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates I=$aln_sorted_BAM O=$aln_sorted_dedup_BAM M=$dedupMetrics_TXT ASSUME_SORTED=true
rm -f $aln_sorted_BAM

### Add read group name
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$aln_sorted_dedup_BAM O=$aln_sorted_dedup_rg_BAM RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=NA12878
rm -f $aln_sorted_dedup_BAM

### Index the BAM file
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar BuildBamIndex I=$aln_sorted_dedup_rg_BAM

###############################################################################
#                      Recalibrate base quality scores                        #
###############################################################################
#set +x;echo "";echo "### Recalibrate base quality scores ###";echo "";set -x

#java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T BaseRecalibrator -R $ref_FASTA -I $aln_sorted_dedup_rg_BAM -knownSites $dbsnps_VCF -knownSites $mills_VCF -knownSites $indels_VCF -o $bqsr_TABLE
#java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T PrintReads -R $ref_FASTA -I $aln_sorted_dedup_rg_BAM -BQSR $bqsr_TABLE -o $aln_sorted_dedup_rg_recal_BAM
#rm -f $aln_sorted_dedup_rg_BAM

###############################################################################
#                      Call variants using Haplotype Caller                   #
###############################################################################
set +x;echo "";echo "### Call variants using Haplotype Caller ###";echo "";set -x

SEC=10
SCC=30
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T HaplotypeCaller -R $ref_FASTA -I $aln_sorted_dedup_rg_BAM --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $raw_variants_VCF

###############################################################################
#                      Recalibrate variant quality scores                     #
###############################################################################
set +x;echo "";echo "### Recalibrate variant quality scores ###";echo "";set -x

resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_VCF"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_VCF"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $1000G_VCF"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_VCF"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
resourceIndels="mills,known=true,training=true,truth=true,prior=12.0 $mills_VCF"
recalParamsIndels="-an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum"
filterLevel="99.0"

java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T VariantRecalibrator -input $raw_variants_VCF -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $snps_RECAL -tranchesFile $snps_TRANCHES -rscriptFile $snps_R
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T ApplyRecalibration -input $raw_variants_VCF -mode SNP -recalFile $snps_RECAL -tranchesFile $snps_TRANCHES --ts_filter_level $filterLevel -o $recalibrated_snps_raw_indels_VCF

java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T VariantRecalibrator -input $recalibrated_snps_raw_indels_VCF -resource:$resourceIndels $recalParamsIndels -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $indels_RECAL -tranchesFile $indels_TRANCHES -rscriptFile $indels_R
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T ApplyRecalibration -input $recalibrated_snps_raw_indels_VCF -mode INDEL -recalFile $indels_RECAL -tranchesFile $indels_TRANCHES --ts_filter_level $filterLevel -o $recalibrated_variants_VCF

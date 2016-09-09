#!/bin/bash

###############################################################################
#     Script for performing GATK Best Practices for Germline SNP & Indel      #
#               Discovery in Whole Genome and Exome Sequence                  #
#     https://software.broadinstitute.org/gatk/best-practices/index.php       #
###############################################################################

### More info: https://software.broadinstitute.org/gatk/events/slides/1506/GATKwr8-A-3-GATK_Best_Practices_and_Broad_pipelines.pdf

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
hapmap_VCF="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_VCF="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
snps_VCF="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_VCF="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_VCF="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_VCF="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

### Programs
install_path="/project/dna/install"
bwa="$install_path/bwa-0.7.13/bwa"
samtools="$install_path/samtools-1.3/bin/samtools"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"

### Temporary directories
javaIOTmpDir="$root/javaIOTmp.dir/"
samtoolsTmpDir="$rootsamtoolTmp.dir/"

### Reference indexing
if [ ! -e "${ref_FASTA}.fai" ]; then
    $samtools faidx $ref_FASTA; date
fi

###############################################################################
#                       Map and mark duplicates                               #
###############################################################################
set +x;echo "";echo "### Map and mark duplicates ###";echo "";set -x

### Generate BWA index
date; $bwa index -a bwtsw $ref_FASTA

### BWA MEM alignment
date; $bwa mem -t $num_threads -M $ref_FASTA $reads_FASTQ > $root.aln_bwa.sam

### Sort SAM file and convert to BAM
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar SortSam I=$root.aln_bwa.sam O=$root.aln_bwa.sorted.bam SORT_ORDER=coordinate
rm -f $root.aln_bwa.sam

### Mark duplicates in the BAM file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates I=$root.aln_bwa.sorted.bam O=$root.aln_bwa.sorted.dupmark.bam M=$root.dedup_metrics.txt ASSUME_SORTED=true
rm -f $root.dedup_metrics.txt
rm -f $root.aln_bwa.sorted.bam

### Add read group name
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$root.aln_bwa.sorted.dupmark.bam O=$root.aln_bwa.sorted.dupmark.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample
rm -f $root.aln_bwa.sorted.dupmark.bam

### Index the BAM file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar BuildBamIndex I=$root.aln_bwa.sorted.dupmark.rg.bam

###############################################################################
#                      Recalibrate base quality scores                        #
###############################################################################
#date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T baseRecalibrator -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.rg.bam -knownSites $dbsnps_VCF -knownSites $mills_VCF -knownSites $indels_VCF -o $root.bqsr.table
#date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T PrintReads -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.rg.bam -BQSR $root.bqsr.table -o $root.aln_bwa.sorted.dupmark.rg.recal.bam
#rm -f $root.bqsr.table
#rm -f $root.aln_bwa.sorted.dupmark.rg.bam

###############################################################################
#                      Call variants using Haplotype Caller                   #
###############################################################################
SEC=10
SCC=30
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T HaplotypeCaller -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.rg.bam --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $root.raw_variants.vcf

###############################################################################
#                      Recalibrate variant quality scores                     #
###############################################################################
resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_VCF"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_VCF"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $1000G_VCF"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_VCF"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
resourceIndels="mills,known=true,training=true,truth=true,prior=12.0 $mills_VCF"
recalParamsIndels="-an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum"
filterLevel="99.0"

date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T VariantRecalibrator -input $root.raw_variants.vcf -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $root.snps.recal -tranchesFile $root.snps.tranches -rscriptFile $root.snps.r
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T ApplyRecalibration -input $root.raw_variants.vcf -mode SNP -recalFile $root.snps.recal -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel -o $root.recalibrated_snps+raw_indels.vcf
rm -f $root.snps.recal
rm -f $root.snps.tranches
rm -f $root.snps.r

date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T VariantRecalibrator -input $root.recalibrated_snps+raw_indels.vcf -resource:$resourceIndels $recalParamsIndels -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $root.indels.recal -tranchesFile $root.indels.tranches -rscriptFile $root.indels.r
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T ApplyRecalibration -input $root.recalibrated_snps+raw_indels.vcf -mode INDEL -recalFile $root.indels.recal -tranchesFile $root.indels.tranches --ts_filter_level $filterLevel -o $root.recalibrated_variants.vcf
rm -f $root.indels.recal
rm -f $root.indels.tranches
rm -f $root.indels.r


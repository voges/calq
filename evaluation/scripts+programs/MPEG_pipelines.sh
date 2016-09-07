#!/bin/bash

###############################################################################
#     Script for performing the MPEG lossy compression framework for          #
#                  genome data (document no. N16324/N100)                     #
###############################################################################

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

### Temporary files
aln_SAM="$base.aln.sam"
aln_BAM="$base.aln.bam"
aln_sorted_BAM="$base.aln.sorted.bam"
aln_sorted_dedup_BAM="$base.aln.sorted.dedup.bam"
aln_sorted_dedup_cleaned_BAM"$base.aln.sorted.dedup.cleaned.bam"
aln_sorted_dedup_cleaned_rg_BAM="$base.aln.sorted.dedup.cleaned.rg.bam"
aln_sorted_dedup_cleaned_rg_realn_BAM="$base.aln.sorted.dedup.cleaned.rg.realn.bam"
aln_sorted_dedup_cleaned_rg_realn_recal_BAM="$base.aln.sorted.dedup.cleaned.rg.realn.recal.bam"
dedupMetrics_TXT="$base.dedupMetrics.txt"
bqsr_TABLE="$base.bqsr.table"
raw_variants_VCF="$base.raw_variants.vcf"
snps_VCF="$base.snps.vcf"
snps_RECAL="$base.snps.recal"
snps_TRANCHES="$base.snps.tranches"
snps_R="$base.snps.R"
filtered_snps_VCF="$base.filtered_snps.vcf"

###############################################################################
#                                 Programs                                    #
###############################################################################
set +x;echo "";echo "### Programs ###";echo "";set -x

install_path="/project/dna/install"
bwa="$install_path/bwa-0.7.13/bwa"
bowtie2="$install_path/bowtie2-2.2.5/bowtie2"
samtools="$install_path/samtools-1.3/bin/samtools"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"

javaIOTmpDir="$base/javaIOTmp.dir"
samtoolsTmpDir="$base/samtoolsTmp.dir"

###############################################################################
#                           Alignment with BWA MEM                            #
###############################################################################
set +x;echo "";echo "### Alignment with BWA MEM ###";echo "";set -x

$bwa index -a bwtsw $ref_FASTA
$bwa mem -t $num_threads -M $ref_FASTA $reads_FASTQ > $aln_bwa_SAM

###############################################################################
#                           Alignment with Bowtie2                            #
###############################################################################
set +x;echo "";echo "### Alignment with Bowtie2 ###";echo "";set -x

$bowtie2-build $ref_FASTA bowtie2_indexes
$bowtie2 -x bowtie2_indexes -1 $reads_FASTQ -S $aln_bowtie2_SAM

###############################################################################
#                             Sorting & indexing                              #
###############################################################################
set +x;echo "";echo "### Sorting & indexing ###";echo "";set -x

### Convert SAM to BAM
$samtools view -@ $num_threads -bh $aln_bwa_SAM > $aln_bwa_BAM
rm -f $aln_bwa_SAM
$samtools view -@ $num_threads -bh $aln_bowtie2_SAM > $aln_bowtie2_BAM
rm -f $aln_bowtie2_SAM

### Sort and index BAM file
$samtools sort -T $samtoolsTmpDir -@ $num_threads -O bam $aln_bwa_BAM > $aln_bwa_sorted_BAM
rm -f $aln_bwa_BAM
$samtools index -T $samtoolsTmpDir -@ $num_threads $aln_bwa_sorted_BAM
$samtools sort -T $samtoolsTmpDir -@ $num_threads -O bam $aln_bowtie2_BAM > $aln_bowtie2_sorted_BAM
rm -f $aln_bowtie2_BAM
$samtools index -T $samtoolsTmpDir -@ $num_threads $aln_bowtie2_sorted_BAM

###############################################################################
#                              Duplicate removal                              #
###############################################################################
set +x;echo "";echo "### Duplicate removal ###";echo "";set -x

### Mark duplicates in the BAM file
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates I=$aln_bwa_sorted_BAM O=$aln_bwa_sorted_dedup_BAM M=$dedupMetrics_bwa_TXT ASSUME_SORTED=true
rm -f $aln_bwa_sorted_BAM
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates I=$aln_bowtie2_sorted_BAM O=$aln_bowtie2_sorted_dedup_BAM M=$dedupMetrics_bowtie2_TXT ASSUME_SORTED=true
rm -f $aln_bowtie2_sorted_BAM

### Remove duplicates
$samtools view -hb -F 0xF40 $aln_bwa_sorted_dedup_BAM > $aln_bwa_sorted_dedup_cleaned_BAM
rm -f $aln_bwa_sorted_dedup_BAM
$samtools view -hb -F 0xF40 $aln_bowtie2_sorted_dedup_BAM > $aln_bowtie2_sorted_dedup_cleaned_BAM
rm -f $aln_bowtie2_sorted_dedup_BAM

### Label the BAM headers and index the resulting file
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$aln_bwa_sorted_dedup_cleaned_BAM O=$aln_bwa_sorted_dedup_cleaned_rg_BAM RGLB=Library RPGL=ILLUMINA RGPU=PlatformUnit RGSM=SampleName
rm -f $aln_bwa_sorted_dedup_cleaned_BAM
$samtools index $aln_bwa_sorted_dedup_cleaned_rg_BAM
java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$aln_bowtie2_sorted_dedup_cleaned_BAM O=$aln_bowtie2_sorted_dedup_cleaned_rg_BAM RGLB=Library RPGL=ILLUMINA RGPU=PlatformUnit RGSM=SampleName
rm -f $aln_bowtie2_sorted_dedup_cleaned_BAM
$samtools index $aln_bowtie2_sorted_dedup_cleaned_rg_BAM

###############################################################################
#                              Indel realignment                              #
###############################################################################
set +x;echo "";echo "### Indel realignment ###";echo "";set -x

java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T RealignerTargetCreator -R $ref_FASTA -I $aln_bwa_sorted_dedup_cleaned_BAM -known $mills_VCF -o $target_INTERVALS
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_FASTA -I $aln_bwa_sorted_dedup_cleaned_BAM -targetIntervals $target_INTERVALS -o $aln_bwa_sorted_dedup_cleaned_rg_realn_BAM
rm -f $aln_bwa_sorted_dedup_cleaned_BAM
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T RealignerTargetCreator -R $ref_FASTA -I $aln_bowtie2_sorted_dedup_cleaned_BAM -known $mills_VCF -o $target_INTERVALS
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_FASTA -I $aln_bowtie2_sorted_dedup_cleaned_BAM -targetIntervals $target_INTERVALS -o $aln_bowtie2_sorted_dedup_cleaned_rg_realn_BAM
rm -f $aln_bowtie2_sorted_dedup_cleaned_BAM


###############################################################################
#                      Base quality score recalibration                       #
###############################################################################
set +x;echo "";echo "### Base quality score recalibration ###";echo "";set -x

java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T BaseRecalibrator -R $ref_FASTA -I $aln_bwa_sorted_dedup_cleaned_rg_realn_BAM -knownSites $dbsnps_VCF -knownSites $mills_VCF -knownSites $indels_VCF -o $bqsr_TABLE
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T PrintReads -R $ref_FASTA -I $aln_bwa_sorted_dedup_cleaned_rg_realn_BAM -BQSR $bqsr_TABLE -o $aln_bwa_sorted_dedup_cleaned_rg_realn_recal_BAM
rm -f $aln_bwa_sorted_dedup_cleaned_rg_realn_BAM
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T BaseRecalibrator -R $ref_FASTA -I $aln_bowtie2_sorted_dedup_cleaned_rg_realn_BAM -knownSites $dbsnps_VCF -knownSites $mills_VCF -knownSites $indels_VCF -o $bqsr_TABLE
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T PrintReads -R $ref_FASTA -I $aln_bowtie2_sorted_dedup_cleaned_rg_realn_BAM -BQSR $bqsr_TABLE -o $aln_bowtie2_sorted_dedup_cleaned_rg_realn_recal_BAM
rm -f $aln_bowtie2_sorted_dedup_cleaned_rg_realn_BAM

###############################################################################
#                          Variant calling with GATK                          #
###############################################################################
set +x;echo "";echo "### Variant calling with GATK ###";echo "";set -x

### Call variants using Haplotype Caller
SEC=10
SCC=30
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T HaplotypeCaller -R $ref_FASTA -I $aln_sorted_dedup_cleaned_rg_realn_recal_BAM --dbnsp $dbsnps_VCF --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $raw_variants_VCF

### SNP extraction
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T SelectVariants -R $ref_FASTA -V $raw_variants_VCF -selectType SNP -o $snps_VCF

### Fiter the variants using VQSR
resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_VCF"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_VCF"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $1000G_VCF"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_VCF"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
filterLevel="99.0"
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T VariantRecalibrator -input $snps_VCF -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $snps_RECAL -tranchesFile $snps_TRANCHES -rscriptFile $snps_R
java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T ApplyRecalibration -input $snps_VCF -mode SNP -recalFile $snps_RECAL -tranchesFile $snps_TRANCHES --ts_filter_level $filterLevel -o $filtered_snps_VCF

###############################################################################
#                        Variant calling with HTSlib                          #
###############################################################################
set +x;echo "";echo "### Variant calling with HTSlib ###";echo "";set -x

### Call variants with SAMtools and BCFtools
#$samtools mpileup -ugf $ref_FASTA $aln_sorted_dedup_cleaned_rg_realn_recal_BAM | $bcftools call -vmO v -o $raw_VCF

### SNP extraction
#java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T SelectVariants -R $ref_FASTA -V $raw_VCF -selectType SNP -o $snps_VCF

### Filter with BCFtools
#$bcftools filter -O v -o $snps_filtered_VCF -s $filtername -i '%QUAL>20' $snps_VCF

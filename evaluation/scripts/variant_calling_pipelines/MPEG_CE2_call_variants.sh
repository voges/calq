#!/bin/bash

###############################################################################
#   Script for performing the ISO/IEC JTC 1/SC 29/WG 11 and ISO/TC 276/WG 5   #
#   Benchmark framework for lossy compression of sequencing quality values    #
#                       (document no. N16727/N156)                            #
#                                    -                                        #
#                             Variant calling                                 #
###############################################################################

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 4 ]; then
    printf "Usage: $0 num_threads log_txt input_bam chromosome\n"
    exit -1
fi

num_threads=$1
printf "Number of threads: $num_threads\n"
log_txt=$2
printf "Log file: $log_txt\n"
input_bam=$3
printf "Input BAM file: $input_bam\n"
chromosome=$4
printf "Chromosome: $chromosome\n"

printf "Checking log file $log_txt ... "
if [ -f $log_txt ]; then printf "log file already exists: $log_txt\n"; exit -1; fi
printf "OK\n"

printf "Checking input BAM file $input_bam ... "
if [ ! -f $input_bam ]; then printf "did not find input BAM file: $input_bam\n"; exit -1; fi
printf "OK\n"

input_bai="$(printf $input_bam | sed 's/\.[^.]*$//')".bai
printf "Checking corresponding BAM index file $input_bai ... "
if [ ! -f $input_bai ]; then printf "did not find BAM index file: $input_bai\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/data/voges/MPEG/GATK_bundle-2.8-b37"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_vcf="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_vcf="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
KG_vcf="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_vcf="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_vcf="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_vcf="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

printf "Checking GATK bundle ... "
if [ ! -f $ref_fasta ]; then printf "did not find $ref_fasta\n"; exit -1; fi
if [ ! -f $hapmap_vcf ]; then printf "did not find $hapmap_vcf\n"; exit -1; fi
if [ ! -f $omni_vcf ]; then printf "did not find $omni_vcf\n"; exit -1; fi
if [ ! -f $KG_vcf ]; then printf "did not find $KG_vcf\n"; exit -1; fi
if [ ! -f $dbsnps_vcf ]; then printf "did not find $dbsnps_vcf\n"; exit -1; fi
if [ ! -f $mills_vcf ]; then printf "did not find $mills_vcf\n"; exit -1; fi
if [ ! -f $indels_vcf ]; then printf "did not find $indels_vcf\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
bowtie2="/project/dna/install/bowtie2-2.2.5/bowtie2"
java="/usr/bin/java"
java_opts=""
samtools="/project/dna/install/samtools-1.3/bin/samtools"

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
picard_jar="/project/dna/install/picard-tools-2.4.1/picard.jar"

printf "Checking executables ... "
if [ ! -x $bowtie2 ]; then printf "did not find $bowtie2\n"; exit -1; fi
if [ ! -x $java ]; then printf "did not find $java\n"; exit -1; fi
if [ ! -x $samtools ]; then printf "did not find $samtools\n"; exit -1; fi
if [ ! -e $GenomeAnalysisTK_jar ]; then printf "did not find $GenomeAnalysisTK_jar\n"; exit -1; fi
if [ ! -e $picard_jar ]; then printf "did not find $picard_jar\n"; exit -1; fi
printf "OK\n"

###############################################################################
#                         Variant calling with GATK                           #
###############################################################################

printf "[1/4] Variant calling ... "
$java $java_opts -jar $GenomeAnalysisTK_jar -T HaplotypeCaller -nct $num_threads -R $ref_fasta -L $chromosome -I $input_bam --dbsnp $dbsnps_vcf --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $input_bam.raw_variants.vcf &>>$log_txt
printf "OK\n"

printf "[2/4] SNP extraction ... "
$java $java_opts -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_fasta -L $chromosome -V $input_bam.raw_variants.vcf -selectType SNP -o $input_bam.snps.vcf &>>$log_txt
printf "OK\n"

printf "[3/4] Variant filtering step 1/2 ... "
resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_vcf"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $KG_vcf"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_vcf"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
$java $java_opts -jar $GenomeAnalysisTK_jar -T VariantRecalibrator -R $ref_fasta -L $chromosome -input $input_bam.snps.vcf -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches -rscriptFile $input_bam.snps.r &>>$log_txt
printf "OK\n"

printf "[4/4] Variant filtering step 2/2 ... "
filterLevel900="90.0"
filterLevel990="99.0"
filterLevel999="99.9"
filterLevel1000="100.0"
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel900 -o $input_bam.snps.filtered900.vcf &>>$log_txt
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel990 -o $input_bam.snps.filtered990.vcf &>>$log_txt
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel999 -o $input_bam.snps.filtered999.vcf &>>$log_txt
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel1000 -o $input_bam.snps.filtered1000.vcf &>>$log_txt
printf "OK\n";

###############################################################################
#                                   Cleanup                                   #
###############################################################################

printf "Cleanup ... "
rm -f $input_bam.snps.r
rm -f $input_bam.snps.recal
rm -f $input_bam.snps.recal.idx
rm -f $input_bam.snps.tranches
printf "OK\n";


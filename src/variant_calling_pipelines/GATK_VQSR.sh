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

if [ "$#" -ne 3 ]; then
    printf "Usage: $0 num_threads input_bam chromosome\n"
    exit -1
fi

num_threads=$1
printf "Number of threads: $num_threads\n"
input_bam=$2
printf "Input BAM file: $input_bam\n"
input_bai="$(printf $input_bam | sed 's/\.[^.]*$//')".bai
printf "Corresponding BAM index file: $input_bai\n"
log_txt="$input.bam.GATK_VQSR.log"
printf "Log file: $log_txt\n"
chromosome=$3
printf "Chromosome: $chromosome\n"

if [ -f $log_txt ]; then printf "Error: File $log_txt file already exists.\n"; exit -1; fi
if [ ! -f $input_bam ]; then printf "Error: Input BAM file $input_bam is not a regular file.\n"; exit -1; fi
if [ ! -f $input_bai ]; then printf "Error: BAM index file $input_bai is not a regular file.\n"; exit -1; fi

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/data/voges/muenteferi/GATK_bundle-2.8-b37"
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
java_opts=""
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
#                         Variant calling with GATK                           #
###############################################################################

printf "[1/4] Variant calling\n"
$java $java_opts -jar $GenomeAnalysisTK_jar -T HaplotypeCaller -nct $num_threads -R $ref_fasta -L $chromosome -I $input_bam --dbsnp $dbsnps_vcf --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $input_bam.GATK.raw_variants.vcf &>>$log_txt

printf "[2/4] SNP extraction\n"
$java $java_opts -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_fasta -L $chromosome -V $input_bam.GATK.raw_variants.vcf -selectType SNP -o $input_bam.GATK.snps.vcf &>>$log_txt

printf "[3/4] Variant filtering step 1/2\n"
resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_vcf"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $KG_vcf"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_vcf"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
$java $java_opts -jar $GenomeAnalysisTK_jar -T VariantRecalibrator -R $ref_fasta -L $chromosome -input $input_bam.GATK.snps.vcf -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches -rscriptFile $input_bam.snps.r &>>$log_txt

printf "[4/4] Variant filtering step 2/2\n"
filterLevel900="90.0"
filterLevel990="99.0"
filterLevel999="99.9"
filterLevel1000="100.0"
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.GATK.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel900 -o $input_bam.GATK.snps.filtered900.vcf &>>$log_txt
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.GATK.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel990 -o $input_bam.GATK.snps.filtered990.vcf &>>$log_txt
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.GATK.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel999 -o $input_bam.GATK.snps.filtered999.vcf &>>$log_txt
$java $java_opts -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_fasta -L $chromosome -input $input_bam.GATK.snps.vcf -mode SNP -recalFile $input_bam.snps.recal -tranchesFile $input_bam.snps.tranches --ts_filter_level $filterLevel1000 -o $input_bam.GATK.snps.filtered1000.vcf &>>$log_txt

###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup\n"
#rm -f $input_bam.snps.r
#rm -f $input_bam.snps.recal
#rm -f $input_bam.snps.recal.idx
#rm -f $input_bam.snps.tranches
printf "Done\n";


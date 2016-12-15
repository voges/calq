#!/bin/bash

###############################################################################
#   Script for performing the ISO/IEC JTC 1/SC 29/WG 11 and ISO/TC 276/WG 5   #
#   Benchmark framework for lossy compression of sequencing quality values    #
#                       (document no. N16525/N119)                            #
#                                    -                                        #
#                             Variant calling                                 #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 num_threads bam_file chromosome"
    exit -1
fi

###############################################################################
#                               Command line                                  #
###############################################################################

date
set -x
script_name=$0
num_threads=$1
bam_file=$2
chromosome=$3

###############################################################################
#                          Data and programs                                  #
###############################################################################

### GATK bundle
gatk_bundle_path="/phys/intern2/tmp/data_gidb/MPEG/GATK_bundle-2.8-b37"
ref_FASTA="$gatk_bundle_path/human_g1k_v37.fasta"
hapmap_VCF="$gatk_bundle_path/hapmap_3.3.b37.vcf"
omni_VCF="$gatk_bundle_path/1000G_omni2.5.b37.vcf"
KG_VCF="$gatk_bundle_path/1000G_phase1.snps.high_confidence.b37.vcf"
dbsnps_VCF="$gatk_bundle_path/dbsnp_138.b37.vcf"
mills_VCF="$gatk_bundle_path/Mills_and_1000G_gold_standard.indels.b37.vcf"
indels_VCF="$gatk_bundle_path/1000G_phase1.indels.b37.vcf"

### Programs
install_path="/project/dna/install"
bowtie2="$install_path/bowtie2-2.2.5/bowtie2"
samtools="$install_path/samtools-1.3/bin/samtools"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"

###############################################################################
#                          Variant calling with GATK                          #
###############################################################################

### Call variants using Haplotype Caller
SEC=10
SCC=30
date; java -jar $GenomeAnalysisTK_jar -T HaplotypeCaller -nct $num_threads -R $ref_FASTA -L $chromosome -I $bam_file --dbsnp $dbsnps_VCF --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $bam_file.raw_variants.vcf

### SNP extraction
date; java -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_FASTA -L $chromosome -V $bam_file.raw_variants.vcf -selectType SNP -o $bam_file.snps.vcf

### Fiter the variants using VQSR
resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_VCF"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_VCF"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $KG_VCF"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_VCF"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
filterLevel900="90.0"
filterLevel990="99.0"
filterLevel999="99.9"
filterLevel1000="100.0"
date; java -jar $GenomeAnalysisTK_jar -T VariantRecalibrator -R $ref_FASTA -L $chromosome -input $bam_file.snps.vcf -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $root.snps.recal -tranchesFile $root.snps.tranches -rscriptFile $root.snps.r
date; java -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_FASTA -L $chromosome -input $bam_file.snps.vcf -mode SNP -recalFile $root.snps.recal -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel900 -o $bam_file.snps.filtered900.vcf
date; java -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_FASTA -L $chromosome -input $bam_file.snps.vcf -mode SNP -recalFile $root.snps.recal -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel990 -o $bam_file.snps.filtered990.vcf
date; java -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_FASTA -L $chromosome -input $bam_file.snps.vcf -mode SNP -recalFile $root.snps.recal -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel999 -o $bam_file.snps.filtered999.vcf
date; java -jar $GenomeAnalysisTK_jar -T ApplyRecalibration -R $ref_FASTA -L $chromosome -input $bam_file.snps.vcf -mode SNP -recalFile $root.snps.recal -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel1000 -o $bam_file.snps.filtered1000.vcf


#!/bin/bash

###############################################################################
#     Script for performing the MPEG lossy compression framework for          #
#                  genome data (document no. N16324/N100)                     #
#                                                                             #
#                            Bowtie2 + GATK HC                                #
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
ref_FASTA="./human_g1k_v37.21.fasta"
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

### Temporary directories
javaIOTmpDir="$root/javaIOTmp.dir/"
samtoolsTmpDir="$root/samtoolsTmp.dir/"

### Reference indexing
if [ ! -e "${ref_FASTA}.fai" ]; then
    date; $samtools faidx $ref_FASTA
fi

###############################################################################
#                           Alignment with Bowtie2                            #
###############################################################################
date; $bowtie2-build $ref_FASTA $root.bowtie2_idx
date; $bowtie2 -x $root.bowtie2_idx -U $reads_FASTQ -S $root.aln_bowtie2.sam --threads $num_threads

###############################################################################
#                             Sorting & indexing                              #
###############################################################################

### Convert SAM to BAM
date; $samtools view -@ $num_threads -bh $root.aln_bowtie2.sam > $root.aln_bowtie2.bam
rm -f $root.aln_bowtie2.sam

### Sort and index BAM file
date; $samtools sort -@ $num_threads -T $samtoolsTmpDir -O bam $root.aln_bowtie2.bam > $root.aln_bowtie2.sorted.bam
rm -f $root.aln_bowtie2.bam
date; $samtools index $root.aln_bowtie2.sorted.bam

###############################################################################
#                              Duplicate removal                              #
###############################################################################

### Mark duplicates in the BAM file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates I=$root.aln_bowtie2.sorted.bam O=$root.aln_bowtie2.sorted.dupmark.bam M=$root.dedup_metrics.txt ASSUME_SORTED=true
rm -f $root.dedup_metrics.txt
rm -f $root.aln_bowtie2.sorted.bam

### Remove duplicates
date; $samtools view -@ $num_threads -bh -F 0xF40 $root.aln_bowtie2.sorted.dupmark.bam > $root.aln_bowtie2.sorted.dupmark.dedup.bam
rm -f $root.aln_bowtie2.sorted.dupmark.bam

### Label the BAM headers and index the resulting file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$root.aln_bowtie2.sorted.dupmark.dedup.bam O=$root.aln_bowtie2.sorted.dupmark.dedup.rg.bam RGID=1 RGLB=Library RGPL=Illumina RGPU=PlatformUnit RGSM=$sample
rm -f $root.aln_bowtie2.sorted.dupmark.dedup.bam
date; $samtools index $root.aln_bowtie2.sorted.dupmark.dedup.rg.bam

###############################################################################
#                              Indel realignment                              #
###############################################################################
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T RealignerTargetCreator -nt $num_threads -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.dedup.rg.bam --known $mills_VCF -o $root.realigner_target.intervals
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.dedup.rg.bam -targetIntervals $root.realigner_target.intervals -o $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.bam
rm -f $root.realigner_target.intervals
rm -f $root.aln_bowtie2.sorted.dupmark.dedup.rg.bam

###############################################################################
#                      Base quality score recalibration                       #
###############################################################################
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T BaseRecalibrator -nct $num_threads -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.bam -knownSites $dbsnps_VCF -knownSites $mills_VCF -knownSites $indels_VCF -o $root.bqsr.table
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T PrintReads -nct $num_threads -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.bam -BQSR $root.bqsr.table -o $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.bam
rm -f $root.bqsr.table
# Keep the BAM file without BQSR!
#rm -f $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.bam

###############################################################################
#                          Variant calling with GATK                          #
###############################################################################

### Call variants using Haplotype Caller
SEC=10
SCC=30
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T HaplotypeCaller -nct $num_threads -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.bam --dbsnp $dbsnps_VCF --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.raw_variants.vcf

### SNP extraction
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T SelectVariants -R $ref_FASTA -V $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.raw_variants.vcf -selectType SNP -o $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.snps.vcf

### Fiter the variants using VQSR
resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_VCF"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_VCF"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $KG_VCF"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_VCF"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
filterLevel="99.0"
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T VariantRecalibrator -input $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.snps.vcf -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $root.snps.recal -tranchesFile $root.snps.tranches -rscriptFile $root.snps.r
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T ApplyRecalibration -input $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.snps.vcf -mode SNP -recalFile $root.snps.recal  -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel -o $root.aln_bowtie2.sorted.dupmark.dedup.rg.realn.recal.snps.filtered.vcf
rm -f $root.snps.recal
rm -f $root.snps.tranches
rm -f $root.snps.r


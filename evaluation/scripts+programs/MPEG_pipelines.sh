#!/bin/bash

###############################################################################
#     Script for performing the MPEG lossy compression framework for          #
#                  genome data (document no. N16324/N100)                     #
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
bowtie2="$install_path/bowtie2-2.2.5/bowtie2"
samtools="$install_path/samtools-1.3/bin/samtools"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"

### Temporary directories
javaIOTmpDir="$root/javaIOTmp.dir"
samtoolsTmpDir="$root/samtoolsTmp.dir"

### Reference indexing
if [ ! -e "${ref_FASTA}.fai" ]; then
    $samtools faidx $ref_FASTA; date
fi

###############################################################################
#                           Alignment with BWA MEM                            #
###############################################################################
date; $bwa index -a bwtsw $ref_FASTA
date; $bwa mem -t $num_threads -M $ref_FASTA $reads_FASTQ > $root.aln_bwa.sam

###############################################################################
#                           Alignment with Bowtie2                            #
###############################################################################
#date; $bowtie2-build $ref_FASTA bowtie2_indexes
#date; $bowtie2 -x bowtie2_indexes -1 $reads_FASTQ -S $root.aln_bowtie2.sam

###############################################################################
#                             Sorting & indexing                              #
###############################################################################

### Convert SAM to BAM
date; $samtools view -@ $num_threads -bh $root.aln_bwa.sam > $root.aln_bwa.bam
rm -f $root.aln_bwa.sam

### Sort and index BAM file
date; $samtools sort -T $samtoolsTmpDir -@ $num_threads -O bam $root.aln_bwa.bam > $root.aln_bwa.sorted.bam
rm -f $root.aln_bwa.bam
date; $samtools index -T $samtoolsTmpDir -@ $num_threads $root.aln_bwa.sorted.bam

###############################################################################
#                              Duplicate removal                              #
###############################################################################

### Mark duplicates in the BAM file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates I=$root.aln_bwa.sorted.bam O=$root.aln_bwa.sorted.dupmark.bam M=$root.dedup_metrics.txt ASSUME_SORTED=true
rm -f $root.aln_bwa.sorted.bam

### Remove duplicates
date; $samtools view -hb -F 0xF40 $root.aln_bwa.sorted.dupmark.bam > $root.aln_bwa.sorted.dupmark.dedup.bam
rm -f $root.aln_bwa.sorted.dupmark.bam

### Label the BAM headers and index the resulting file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$root.aln_bwa.sorted.dupmark.dedup.bam O=$root.aln_bwa.sorted.dupmark.dedup.rg.bam RGLB=Library RPGL=ILLUMINA RGPU=PlatformUnit RGSM=SampleName
rm -f $root.aln_bwa.sorted.dupmark.dedup.bam
date; $samtools index $root.aln_bwa.sorted.dupmark.dedup.rg.bam

###############################################################################
#                              Indel realignment                              #
###############################################################################
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T RealignerTargetCreator -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.dedup.rg.bam -known $mills_VCF -o $root.realigner_target.intervals
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.dedup.rg.bam -targetIntervals $root.realigner_target.intervals -o $root.aln_bwa.sorted.dupmark.dedup.rg.realn.bam
rm -f $root.realigner_target.intervals
rm -f $root.aln_bwa.sorted.dupmark.dedup.rg.bam

###############################################################################
#                      Base quality score recalibration                       #
###############################################################################
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T baseRecalibrator -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.dedup.rg.realn.bam -knownSites $dbsnps_VCF -knownSites $mills_VCF -knownSites $indels_VCF -o $root.bqsr.table
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -nct $num_threads -T PrintReads -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.dedup.rg.realn.bam -BQSR $root.bqsr.table -o $root.aln_bwa.sorted.dupmark.dedup.rg.realn.recal.bam
rm -f $root.bqsr.table
rm -f $root.aln_bwa.sorted.dupmark.dedup.rg.realn.bam

###############################################################################
#                          Variant calling with GATK                          #
###############################################################################

### Call variants using Haplotype Caller
SEC=10
SCC=30
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T HaplotypeCaller -R $ref_FASTA -I $root.aln_bwa.sorted.dupmark.dedup.rg.realn.recal.bam --dbnsp $dbsnps_VCF --genotyping_mode DISCOVERY -stand_emit_conf $SEC -stand_call_conf $SCC -o $root.raw_variants.vcf

### SNP extraction
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T SelectVariants -R $ref_FASTA -V $root.raw_variants.vcf -selectType SNP -o $root.snps.vcf

### Fiter the variants using VQSR
resourceSNPs1="hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_VCF"
resourceSNPs2="omni,known=false,training=true,truth=true,prior=12.0 $omni_VCF"
resourceSNPs3="1000G,known=false,training=true,truth=false,prior=10.0 $1000G_VCF"
resourceSNPs4="dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnps_VCF"
recalParamsSNPs="-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum"
filterLevel="99.0"
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T VariantRecalibrator -input $root.snps.vcf -resource:$resourceSNPs1 -resource:$resourceSNPs2 -resource:$resourceSNPs3 -resource:$resourceSNPs4 $recalParamsSNPs -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $root.snps.recal -tranchesFile $root.snps.tranches -rscriptFile $root.snps.r
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -R $ref_FASTA -T ApplyRecalibration -input $root.snps.vcf -mode SNP -recalFile $root.snps.recal  -tranchesFile $root.snps.tranches --ts_filter_level $filterLevel -o $root.snps.filtered.vcf
rm -f $root.snps.recal
rm -f $root.snps.tranches
rm -f $root.snps.r

###############################################################################
#                        Variant calling with HTSlib                          #
###############################################################################

### Call variants
#date; $samtools mpileup -ugf $ref_FASTA $root.aln_bwa.sorted.dupmark.dedup.rg.realn.recal.bam | $bcftools call -vmO v -o $root.raw_variants.vcf

### SNP extraction
#date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T SelectVariants -R $ref_FASTA -V $root.raw_variants.vcf -selectType SNP -o $root.snps.vcf

### Filter with BCFtools
#date; $bcftools filter -O v -o $root.snps.filtered.vcf -s $filtername -i '%QUAL>20' $root.snps.vcf

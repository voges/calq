#!/bin/bash

###############################################################################
#   Script for performing the ISO/IEC JTC 1/SC 29/WG 11 and ISO/TC 276/WG 5   #
#   Benchmark framework for lossy compression of sequencing quality values    #
#                       (document no. N16525/N119)                            #
#                                    -                                        #
#                                Alignment                                    #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 num_threads reads pairing platform"
    exit -1
fi

###############################################################################
#                               Command line                                  #
###############################################################################

date
set -x
script_name=$0
num_threads=$1
reads=$2
pairing=$3 # paired or unpaired
platform=$4 # ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
root="$reads"

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

### Temporary directories
javaIOTmpDir="$root/javaIOTmp.dir/"

###############################################################################
#                           Alignment with Bowtie2                            #
###############################################################################

date; $bowtie2-build $ref_FASTA $root.bowtie2_idx
if [ "$pairing" = "unpaired" ]; then
    date; $bowtie2 -x $root.bowtie2_idx -U $reads.fastq -S $root.aln_bowtie2.sam --threads $num_threads
else
    if [ "$pairing" = "paired" ]; then
        date; $bowtie2 -x $root.bowtie2_idx -1 $reads\_1.fastq -2 $reads\_2.fastq -S $root.aln_bowtie2.sam --threads $num_threads
    else
        echo "pairing argument must be either 'unpaired' or 'paired'"
        exit -1
    fi
fi
rm -f $root.bowtie2_idx*

###############################################################################
#                             Sorting & indexing                              #
###############################################################################

### Convert SAM to BAM
date; $samtools view -@ $num_threads -bh $root.aln_bowtie2.sam > $root.aln_bowtie2.bam
rm -f $root.aln_bowtie2.sam

### Sort and index BAM file
date; $samtools sort -@ $num_threads -O bam $root.aln_bowtie2.bam > $root.aln_bowtie2.sorted.bam
rm -f $root.aln_bowtie2.bam
date; $samtools index $root.aln_bowtie2.sorted.bam

###############################################################################
#                              Duplicate removal                              #
###############################################################################

### Mark duplicates in the BAM file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates I=$root.aln_bowtie2.sorted.bam O=$root.aln_bowtie2.sorted.dupmark.bam M=$root.dedup_metrics.txt ASSUME_SORTED=true
rm -f $root.dedup_metrics.txt
rm -f $root.aln_bowtie2.sorted.bam
rm -f $root.aln_bowtie2.sorted.bam.bai

### Label the BAM headers and index the resulting file
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar AddOrReplaceReadGroups I=$root.aln_bowtie2.sorted.dupmark.bam O=$root.aln_bowtie2.sorted.dupmark.rg.bam RGID=1 RGLB=Library RGPL=$platform RGPU=PlatformUnit RGSM=$root
rm -f $root.aln_bowtie2.sorted.dupmark.bam
date; $samtools index $root.aln_bowtie2.sorted.dupmark.rg.bam

###############################################################################
#                              Indel realignment                              #
###############################################################################

date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T RealignerTargetCreator -nt $num_threads -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.rg.bam --known $mills_VCF -o $root.realigner_target.intervals
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.rg.bam -targetIntervals $root.realigner_target.intervals -o $root.aln_bowtie2.sorted.dupmark.rg.realn.bam
rm -f $root.realigner_target.intervals
rm -f $root.aln_bowtie2.sorted.dupmark.rg.bam
rm -f $root.aln_bowtie2.sorted.dupmark.rg.bam.bai

###############################################################################
#                      Base quality score recalibration                       #
###############################################################################

date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T BaseRecalibrator -nct $num_threads -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.rg.realn.bam -knownSites $dbsnps_VCF -knownSites $mills_VCF -knownSites $indels_VCF -o $root.bqsr.table
date; java -jar -Djava.io.tmpdir=$javaIOTmpDir $GenomeAnalysisTK_jar -T PrintReads -nct $num_threads -R $ref_FASTA -I $root.aln_bowtie2.sorted.dupmark.rg.realn.bam -BQSR $root.bqsr.table -o $root.aln_bowtie2.sorted.dupmark.rg.realn.recal.bam
rm -f $root.bqsr.table
rm -f $root.aln_bowtie2.sorted.dupmark.rg.realn.bam
rm -f $root.aln_bowtie2.sorted.dupmark.rg.realn.bai

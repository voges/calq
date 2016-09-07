#!/bin/bash

###############################################################################
#    Script for performing WGS/WES Mapping to Variant Calls - Version 1.0     #
#                      http://www.htslib.org/workflow/                        #
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

### File names
aln_SAM=$base".aln.sam"
aln_fixmate_BAM=$base".aln.fixmate.bam"
aln_fixmate_sorted_BAM=$base".aln.fixmate.sorted.bam"
aln_fixmate_sorted_INTERVALS=$base".aln.fixmate.sorted.intervals"
aln_fixmate_sorted_realn_BAM=$base".aln.fixmate.sorted.realn.bam"
aln_fixmate_sorted_realn_dedup_BAM=$base".aln.fixmate.sorted.realn.dedup.bam"
raw_variants_VCF_GZ=$base".raw_variants.vcf.gz"
filtered_variants_VCF_GZ=$base".filtered_variants.vcf.gz"

###############################################################################
#                                Programs                                     #
###############################################################################
set +x;echo "";echo "### Programs ###";echo "";set -x

install_path="/project/dna/install"
bwa="$install_path/bwa-0.7.13/bwa"
samtools="$install_path/samtools-1.3/bin/samtools"
picard_jar="$install_path/picard-tools-2.4.1/picard.jar"
GenomeAnalysisTK_jar="$install_path/gatk-3.6/GenomeAnalysisTK.jar"
tabix="$install_path/project/dna/install/htslib-1.3/bin/tabix"

### Temporary directories
javaIOTmpDir="$base/javaIOTmp.dir"
samtoolsTmpDir="$base/samtoolsTmp.dir"

###############################################################################
#                                  Mapping                                    #
###############################################################################
set +x;echo "";echo "### Mapping ###";echo "";set -x

### Generate BWA index
$bwa index -a bwtsw $ref_FASTA

### BWA MEM alignment
$bwa mem -t $num_threads -R '@RG\tID:1\tSM:NA12878\tLB:lib1' -M $ref_FASTA $reads_FASTQ > $aln_SAM

### Clean up read pairing information and flags and convert to BAM
$samtools fixmate -O bam $aln_SAM $aln_fixmate_BAM
rm -f $aln_SAM

### Sort BAM file
$samtools sort -T $samtoolsTmpDir -@ $num_threads -O bam $aln_fixmate_BAM > $aln_fixmate_sorted_BAM
rm -f $aln_fixmate_BAM

###############################################################################
#                                 Improvement                                 #
###############################################################################
set +x;echo "";echo "### Improvement ###";echo "";set -x

### Reduce the number of miscalls of INDELs by realigning
java -Xmx2g -jar $GenomeAnalysisTK_jar -T RealignerTargetCreator -R $ref_FASTA -I $aln_fixmate_sorted_BAM -o $aln_fixmate_sorted_INTERVALS --known $mills_VCF

java -Xmx4g -jar $GenomeAnalysisTK_jar -T IndelRealigner -R $ref_FASTA -I $aln_fixmate_sorted_BAM -targetIntervals $aln_fixmate_sorted_INTERVALS --known $mills_VCF -o $aln_fixmate_sorted_realn_BAM
rm -f $aln_fixmate_sorted_BAM

### Mark duplicates in the BAM file
java -Xmx2g -jar -Djava.io.tmpdir=$javaIOTmpDir $picard_jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=$aln_fixmate_sorted_realn_BAM O=$aln_fixmate_sorted_realn_dedup_BAM
rm -f $aln_fixmate_sorted_realn_BAM

### Index the BAM file
$samtools index $aln_fixmate_sorted_realn_dedup_BAM

###############################################################################
#                               Variant calling                               #
###############################################################################
set +x;echo "";echo "### Variant calling ###";echo "";set -x

### Call variants
$samtools mpileup -ugf $ref_FASTA $aln_fixmate_sorted_realn_dedup_BAM | $bcftools call -vmO z -o $raw_variants_VCF_GZ

### Index VCF file
$tabix -p vcf $raw_variants_VCF_GZ

### Filter variants
$bcftools filter -O v -o $filtered_variants_VCF_GZ -s LOWQUAL -i '%QUAL>20' $raw_variants_VCF_GZ

#!/bin/bash

###############################################################################
#      Script for performing SNP calling and hard filtering with DeepVariant  #
#                                                                             #
#   More info at http://gatkforums.broadinstitute.org/gatk/discussion/2806/   #
###############################################################################

###############################################################################
#                               Command line                                  #
###############################################################################

if [ "$#" -ne 4 ]; then
    printf "Usage: $0 num_threads input_bam chromosome GPU_ID\n"
    exit -1
fi

num_threads=$1
printf "Number of threads: $num_threads\n"
input_bam=$2
printf "Input BAM file: $input_bam\n"
input_bai="$(printf $input_bam | sed 's/\.[^.]*$//')".bai
printf "Corresponding BAM index file: $input_bai\n"
log_txt="$input.bam.DeepVariant.log"
printf "Log file: $log_txt\n"
chromosome=$3
printf "Chromosome: $chromosome\n"
gpu=$4
printf "GPU: $gpu\n"

#Select GPU
export CUDA_VISIBLE_DEVICES=$gpu

if [ -f $log_txt ]; then printf "Error: File $log_txt file already exists.\n"; exit -1; fi
if [ ! -f $input_bam ]; then printf "Error: Input BAM file $input_bam is not a regular file.\n"; exit -1; fi
if [ ! -f $input_bai ]; then printf "Error: BAM index file $input_bai is not a regular file.\n"; exit -1; fi

###############################################################################
#                                GATK bundle                                  #
###############################################################################

gatk_bundle_path="/data/voges/muenteferi/GATK_bundle-2.8-b37"
dv_bundle="/home/muenteferi/Dokumente/deepvariant"
dv_examples="$dv_bundle/bin/make_examples.zip"
dv_calling="$dv_bundle/bin/call_variants.zip"
dv_postprocess="$dv_bundle/bin/postprocess_variants.zip"
ref_fasta="$gatk_bundle_path/human_g1k_v37.fasta"
dv_model="$dv_bundle/model/model.ckpt"

if [ ! -f $ref_fasta ]; then printf "Error: File $ref_fasta is not a regular file.\n"; exit -1; fi

###############################################################################
#                                Executables                                  #
###############################################################################

# Binaries
java="/usr/bin/java"
java_opts=""

# JAR files
GenomeAnalysisTK_jar="/project/dna/install/gatk-3.6/GenomeAnalysisTK.jar"
python="/home/muenteferi/Dokumente/deepvariant/examplepy/bin/python"

if [ ! -x $java ]; then printf "Error: Binary file $java is not executable.\n"; exit -1; fi
if [ ! -f $GenomeAnalysisTK_jar ]; then printf "Error: JAR file $GenomeAnalysisTK_jar is not a regular file.\n"; exit -1; fi
if [ ! -f $dv_examples ]; then printf "Error: JAR file $dv_examples is not a regular file.\n"; exit -1; fi
if [ ! -f $dv_calling ]; then printf "Error: JAR file $dv_calling is not a regular file.\n"; exit -1; fi
if [ ! -f $dv_postprocess ]; then printf "Error: JAR file $dv_postprocess is not a regular file.\n"; exit -1; fi

###############################################################################
#                         Variant calling with GATK                           #
###############################################################################

LOGDIR=./logs
N_SHARDS=$num_threads

mkdir -p "${LOGDIR}"

printf "[1/4] Extract candidate sites\n"

if [ ! -f $input_bam.examples.tfrecord-00000-of-0000$N_SHARDS.gz ] 
then

time seq 0 $((N_SHARDS-1)) | parallel --eta --halt 2 --joblog "${LOGDIR}/log" --res "${LOGDIR}" \
  $python $dv_examples \
    --mode calling \
    --ref "$ref_fasta" \
    --reads "$input_bam" \
    --examples "$input_bam.examples.tfrecord@${N_SHARDS}.gz" \
    --regions "$chromosome" \
    --task {}
else
printf "Examples $input_bam.examples existing, skipping make_examples.\n"
fi
printf "[2/4] Variant calling\n"
echo "$python $dv_calling --outfile \"$input_bam.examples.tfrecord_filtered.gz\"  --examples \"$input_bam.examples.tfrecord@${N_SHARDS}.gz\"  --checkpoint \"$dv_model\" &>>$log_txt"
$python $dv_calling --outfile "$input_bam.examples.tfrecord_filtered.gz"  --examples "$input_bam.examples.tfrecord@${N_SHARDS}.gz"  --checkpoint "$dv_model" &>>$log_txt

printf "[3/4] Postprocessing\n"
$python $dv_postprocess --ref $ref_fasta --infile "$input_bam.examples.tfrecord_filtered.gz" --outfile "$input_bam.DeepVariant.raw_variants.vcf" &>>$log_txt

printf "[4/4] SNP / Indel extraction\n"
$java $java_opts -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_fasta -L $chromosome -V "$input_bam.DeepVariant.raw_variants.vcf" -selectType SNP -o "$input_bam.DeepVariant.snps.vcf" &>>$log_txt
$java $java_opts -jar $GenomeAnalysisTK_jar -T SelectVariants -R $ref_fasta -L $chromosome -V "$input_bam.DeepVariant.raw_variants.vcf" -selectType INDEL -o "$input_bam.DeepVariant.indels.vcf" &>>$log_txt



###############################################################################
#                                   Cleanup                                   #
###############################################################################

#printf "Cleanup ... "
#
printf "Done\n";


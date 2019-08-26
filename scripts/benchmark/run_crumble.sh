#!/usr/bin/env bash

# Run Crumble's BAM-to-CRAM compression

readonly cram_size=""
readonly crumble=""
readonly samtools=""
readonly scramble=""
readonly time="$(whereis time)"

if [[ "${#}" -ne 3 ]]; then
    echo "Usage: ${0} input_bam num_threads"
    exit 1
fi

readonly input_bam="${1}"
readonly ref_fasta="${2}"
readonly num_threads="${3}"

echo "Constructing FASTA index file: ${ref_fasta}.fai"
if [[ -f "${ref_fasta}.fai" ]]; then
    echo "${ref_fasta}.fai already exists (not reproducing it)"
else
    "${samtools}" faidx "${ref_fasta}"
fi

printf "Crumble BAM-to-BAM encoding with compression level 1\n"
cmd="${crumble} -v -1 ${input_bam} ${input_bam}.${crumble_string}-1.bam"
${time} -v -o ${input_bam}.$crumble_string-1.time ${cmd} &> ${input_bam}.${crumble_string}-1.log

printf "BAM index creation\n"
${samtools} index -b ${input_bam}.$crumble_string-1.bam ${input_bam}.${crumble_string}-1.bai

printf "Scramble BAM-to-CRAM encoding\n"
${scramble} -r ${ref_fasta} -t ${num_threads} ${input_bam}.${crumble_string}-1.bam ${input_bam}.${crumble_string}-1.bam.cram &> ${input_bam}.${crumble_string}-1.bam.${scramble_string}.log

printf "Reporting CRAM size\n"
${cram_size} ${input_bam}.${crumble_string}-1.bam.cram &> ${input_bam}.${crumble_string}-1.bam.${scramble_string}.cram_size
mv ${input_bam}.${crumble_string}-1.bam.cram ${input_bam}.${crumble_string}-1.bam.${scramble_string}

printf "Crumble BAM-to-BAM encoding with compression level 9\n"
cmd="${crumble} -v -9 ${input_bam} ${input_bam}.${crumble_string}-9.bam"
${time} -v -o ${input_bam}.${crumble_string}-9.time ${cmd} &> ${input_bam}.${crumble_string}-9.log

printf "BAM index creation\n"
${samtools} index -b ${input_bam}.${crumble_string}-9.bam ${input_bam}.${crumble_string}-9.bai

printf "Scramble BAM-to-CRAM encoding\n"
${scramble} -r ${ref_fasta} -t ${num_threads} ${input_bam}.${crumble_string}-9.bam ${input_bam}.${crumble_string}-9.bam.cram &> ${input_bam}.${crumble_string}-9.bam.${scramble_string}.log

printf "Reporting CRAM size\n"
${cram_size} ${input_bam}.${crumble_string}-9.bam.cram &> ${input_bam}.${crumble_string}-9.bam.${scramble_string}.cram_size
mv ${input_bam}.${crumble_string}-9.bam.cram ${input_bam}.${crumble_string}-9.bam.${scramble_string}

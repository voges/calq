#!/usr/bin/env bash

# Extract the quality values from a FASTQ file and compress them with bgzip

readonly bgzip=""
readonly extract_part_fastq_py=""
readonly python=""

if [[ "${#}" -ne 2 ]]; then
    echo "Usage: ${0} input_fastq num_threads"
    exit 1
fi

readonly input_fastq="${1}"
readonly num_threads="${2}"

"${python}" "${extract_part_fastq_py}" "${input_fastq}" 3 1>"${input_fastq}.qual"

"${bgzip}" -@ "${num_threads}" -c "${input_fastq}.qual" 1>"${input_fastq}.qual.bgz"

#!/usr/bin/env bash

# Extract a chromosome from a BED file

if [[ "${#}" -ne 2 ]]; then
    echo "Usage: ${0} <file.bed> <chromosome>"
    exit 1
fi
readonly output_prefix=${1%.*} # strip .bed
readonly chromosome=${2}
if [[ ! -f "${1}" ]]; then
    echo "Error: '${1}' is not a regular file"
    exit 1
fi

readonly output_file="${output_prefix}.${chromosome}.bed"
if [[ -f "${output_file}" ]]; then
    echo "'${output_file}' already exists (not reproducing it)"
    exit 0
fi
echo "Writing output to: ${output_file}"
grep -P "^${chromosome}" "${1}" 1>"${output_file}"

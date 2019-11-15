#!/usr/bin/env bash

# Extract a chromosome from a SAM or BAM file

# Required programs
readonly samtools="/project/omics/install/samtools-1.3/bin/samtools"

# Command line
if [[ "${#}" -ne 2 ]]; then
    echo "Usage: ${0} <file.[sam|bam]> <chromosome>"
    exit 1
fi
readonly input_sambam=${1}
readonly output_prefix=${1%.*} # strip .sam/.bam
readonly chromosome=${2}
if [[ ! -f "${input_sambam}" ]]; then
    echo "Error: '${input_sambam}' is not a regular file"
    exit 1
fi

# Convert input file from SAM to BAM, if necessary
if [[ ${1} == *.sam ]]; then
    echo "SAM-to-BAM conversion: ${input_sambam} -> ${output_prefix}.bam"
    if [[ -f "${output_prefix}.bam" ]]; then
        echo "BAM file '${output_prefix}.bam' already exists (not reproducing it)"
    else
        ${samtools} view -bh "${1}" 1>"${output_prefix}.bam" || exit 1
    fi
fi

# Construct BAM index file, if necessary
echo "Constructing BAM index file: ${output_prefix}.bai"
if [[ -f "${output_prefix}.bai" ]]; then
    echo "BAM index file '${output_prefix}.bai' already exists (not reproducing it)"
else
    ${samtools} index "${output_prefix}.bam" "${output_prefix}.bai" || exit 1
fi

# Extract the chromosome
echo "Extracting chromosome: ${chromosome}"
if [[ -f "${output_prefix}.${chromosome}.sam" ]]; then
    echo "SAM file '${output_prefix}.${chromosome}.sam' already exists (not reproducing it)"
else
    ${samtools} view -h "${output_prefix}.bam" "${chromosome}" 1>"${output_prefix}.${chromosome}.sam" || exit 1
fi

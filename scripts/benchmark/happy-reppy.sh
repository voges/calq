#!/usr/bin/env bash

# Script to report variant calling results with
#   hap.py (https://github.com/Illumina/hap.py)
# and
#   rep.py (https://github.com/ga4gh/benchmarking-tools/)

set -euo pipefail

readonly python=""
readonly hap_py=""
readonly rep_py=""

if [[ "${#}" -ne 6 ]]; then
    echo "Usage: ${0} num_threads log_txt variants_vcf ref_fasta golden_vcf_gz golden_bed"
    exit 1
fi

readonly num_threads="${1}"
readonly log_txt="${2}"
readonly variants_vcf="${3}"
readonly ref_fasta="${4}"
readonly golden_vcf_gz="${5}"
readonly golden_bed="${6}"

echo "[1/2] Running hap.py"
"${python}" "${hap_py}" \
    --verbose \
    --threads "${num_threads}" \
    "${golden_vcf_gz}" \
    "${variants_vcf}" \
    -f ${golden_bed} \
    -o ${variants_vcf}.happy \
    -r ${ref_fasta} \
    --roc VQLSOD \
    &>>${log_txt}

echo "[2/2] Running rep.py"
rm -f "${variants_vcf}.happy.rep.tsv"
echo "method\tcomparisonmethod\tfiles" >>"${variants_vcf}".happy.rep.tsv
echo "${variants_vcf}\t${golden_vcf_gz}\t" >>"${variants_vcf}".happy.rep.tsv
echo "${variants_vcf}.happy.extended.csv," >>"${variants_vcf}".happy.rep.tsv
for i in "${variants_vcf}".happy.roc.Locations.*; do
    echo "${i},";
done >> "${variants_vcf}".happy.rep.tsv
sed -i '$ s/.$//' "${variants_vcf}".happy.rep.tsv
"${python}" "${rep_py}" \
    -o "${variants_vcf}.happy.rep.tsv.reppy.html" \
    -l "${variants_vcf}.happy.rep.tsv" \
    &>>${log_txt}

#!/usr/bin/env bash

set -euxo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

readonly build_dir="${git_root_dir}/cmake-build-debug-all"
if [[ ! -d "${build_dir}" ]]; then
    echo "Error: build directory does not exist: ${build_dir}"
    exit 1
fi

# Check wether any .gcda files are present
num_gcda_files=$(find "${build_dir}" -type f -name "*.gcda" | wc -l)
num_gcda_files="${num_gcda_files//[[:space:]]/}"
if [[ "${num_gcda_files}" == 0 ]]; then
    echo "Error: no .gcda files found in build directory: ${build_dir}"
    exit 1
fi

cd "${git_root_dir}"

# Capture coverage info
lcov --directory . --capture --output-file coverage.info

# Filter out files
lcov --remove coverage.info \
    --output-file coverage.info \
    '/usr/*' \
    '*/calq/cmake-build-*/*' \
    '*/calq/tests/*'

# Output coverage data on the console (optional)
lcov --list coverage.info

# Upload report to codevio.io (if run on Travis) or generate local HTML report
set +u # in the following lines do *not* treat unset variables as an error
# shellcheck disable=SC2236
if [[ ! -z "${CI}" ]]; then
    bash <(curl -s https://codecov.io/bash) -f coverage.info -t "${CODECOV_TOKEN}"
else
    readonly local_codecov_dir="${build_dir}/codecov/html"
    genhtml coverage.info --output-directory "${local_codecov_dir}"
    set +x; echo ""; echo "Coverage report generated locally at: ${local_codecov_dir}"; echo ""; set -x
    set +x; echo ""; echo "Maybe you'd like to use this command: firefox ${local_codecov_dir}/index.html &"; echo ""; set -x

fi
set -u # treat unset variables as error again

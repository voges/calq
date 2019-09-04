#!/usr/bin/env bash

set -euxo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

readonly build_dir="${git_root_dir}/cmake-build-debug-with-doc-and-coverage"
[[ -d "${build_dir}" ]] # exit if build dir does *not* exist

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

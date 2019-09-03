#!/usr/bin/env bash

set -euxo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

readonly build_dir="${git_root_dir}/cmake-build-debug-with-doc"
[[ -d "${build_dir}" ]] # exit if build dir does *not* exist

cd "${git_root_dir}"

# Capture coverage info
lcov --directory . --capture --output-file coverage.info

# Filter out files
lcov --remove coverage.info \
    '/usr/*' \
    '*/calq/build/*' \
    '*/calq/cmake-build/*' \
    '*/calq/cmake-build-debug/*' \
    '*/calq/cmake-build-debug-with-doc/*' \
    '*/calq/cmake-build-release/*' \
    '*/calq/tests/*' \
    --output-file coverage.info

# Output coverage data on the console (optional)
lcov --list coverage.info

# Upload report to codevio.io (if run on Travis) or generate local HTML report
set +u # in the following lines do *not* treat unset variables as an error
if [[ ! -z "${CI}" ]]; then
    bash <(curl -s https://codecov.io/bash) -f coverage.info -t "${CODECOV_TOKEN}"
else
    genhtml coverage.info --output-directory "${build_dir}/codecov/html/"
fi
set -u # treat unset variables as error again

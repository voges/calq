#!/usr/bin/env bash

set -euxo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

cd "${git_root_dir}"

# Capture coverage info
lcov --directory . --capture --output-file coverage.info

# Filter out files
lcov --remove coverage.info '*/calq/build/*' '/usr/*' '*/calq/tests/*' --output-file coverage.info

# Output coverage data on the console (optional)
lcov --list coverage.info

# Generate HTML output
#genhtml coverage.info --output-directory ./build/codecov/html/

# Upload report to codevio.io
bash <(curl -s https://codecov.io/bash) -f coverage.info -t "${CODECOV_TOKEN}" || echo "Codecov did not collect coverage reports"

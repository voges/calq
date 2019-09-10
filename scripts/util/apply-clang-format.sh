#!/usr/bin/env bash

# Apply clang-format in-place on all C++ source code files

set -euo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

dirs=()
dirs+=("${git_root_dir}/source/apps/cip")
dirs+=("${git_root_dir}/source/libs/calq")
dirs+=("${git_root_dir}/tests/example")
dirs+=("${git_root_dir}/tests/libs/calq")

extensions=()
extensions+=("h")
extensions+=("cc")

cd "${git_root_dir}"

for i in "${!dirs[@]}"; do
    source_dir=${dirs[${i}]}
    echo "Processing directory: ${source_dir}"
    for j in "${!extensions[@]}"; do
        extension=${extensions[${j}]}
        echo "Processing all *.${extension} files"
        for file in "${source_dir}"/*."${extension}"; do
            if [[ -f "${file}" ]]; then
                echo "Running clang-format on: ${file}"
                clang-format -i "${file}"
            else
                echo "Skipping: ${file}"
            fi
        done
    done
done

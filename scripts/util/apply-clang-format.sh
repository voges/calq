#!/usr/bin/env bash

# Apply clang-format in-place on all C++ source code files

set -euo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

cd "${git_root_dir}"

for f in "${git_root_dir}"/source/apps/cip/*.h; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

for f in "${git_root_dir}"/source/apps/cip/*.cc; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

for f in "${git_root_dir}"/source/libs/calq/*.h; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

for f in "${git_root_dir}"/source/libs/calq/*.cc; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

for f in "${git_root_dir}"/source/libs/util/*.h; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

for f in "${git_root_dir}"/source/libs/util/*.cc; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

# for f in "${git_root_dir}"/tests/example/*.h; do
#     echo "Running clang-format on: ${f}"
#     clang-format -i "${f}"
# done

for f in "${git_root_dir}"/tests/example/*.cc; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

# for f in "${git_root_dir}"/tests/libs/calq/*.h; do
#     echo "Running clang-format on: ${f}"
#     clang-format -i "${f}"
# done

for f in "${git_root_dir}"/tests/libs/calq/*.cc; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

# for f in "${git_root_dir}"/tests/libs/util/*.h; do
#     echo "Running clang-format on: ${f}"
#     clang-format -i "${f}"
# done

for f in "${git_root_dir}"/tests/libs/util/*.cc; do
    echo "Running clang-format on: ${f}"
    clang-format -i "${f}"
done

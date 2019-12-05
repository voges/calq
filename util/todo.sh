#!/usr/bin/env bash

set -euo pipefail

self="${0}"
self_name="${self##*/}"

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

dirs=()
dirs+=("${git_root_dir}/src")
dirs+=("${git_root_dir}/test")

for i in "${!dirs[@]}"; do
    dir=${dirs[${i}]}
    echo "[${self_name}] listing TODOs and FIXMEs in tree: ${dir}"
    set +e
    grep --color --recursive --line-number --word-regexp 'TODO' "${dir}"
    grep --color --recursive --line-number --word-regexp 'FIXME' "${dir}"
    set -e
done

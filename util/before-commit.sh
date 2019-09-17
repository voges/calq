#!/usr/bin/env bash

set -euo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

execute() {
    echo ""
    echo "==============================================================================="
    echo "Running: ${1}"
    echo "==============================================================================="
    echo ""

    "${1}"
}

execute "${git_root_dir}/util/clang-format.sh"
execute "${git_root_dir}/util/authors.sh"
execute "${git_root_dir}/util/shellcheck.sh"

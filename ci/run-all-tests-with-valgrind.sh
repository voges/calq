#!/usr/bin/env bash

set -euxo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

readonly build_dir="${git_root_dir}/cmake-build-debug-all"
readonly calq_tests_app="${build_dir}/bin/calq-tests"
[[ -x "${calq_tests_app}" ]] # exit if calq_tests_app does *not* exist

valgrind \
    --tool=memcheck \
    --verbose \
    --track-origins=yes \
    --leak-check=full \
    --show-leak-kinds=all \
    "${calq_tests_app}"

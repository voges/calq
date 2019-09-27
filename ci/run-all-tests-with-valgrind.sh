#!/usr/bin/env bash

set -euxo pipefail

self="${0}"
self_name="${self##*/}"

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

readonly build_dir="${git_root_dir}/cmake-build-debug-all"
readonly calq_test_app="${build_dir}/bin/calq-test"
if [[ ! -x "${calq_test_app}" ]]; then
    echo "[${self_name}] error: calq-test application does not exist: ${calq_test_app}"
    exit 1
fi

valgrind \
    --tool=memcheck \
    --verbose \
    --track-origins=yes \
    --leak-check=full \
    --show-leak-kinds=all \
    "${calq_test_app}"

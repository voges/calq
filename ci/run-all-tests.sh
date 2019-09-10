#!/usr/bin/env bash

set -euxo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

readonly build_dir="${git_root_dir}/cmake-build-debug-with-doc-and-coverage"
[[ -d "${build_dir}" ]] # exit if build dir does *not* exist

"${build_dir}/bin/calq-tests"

#!/usr/bin/env bash

set -euxo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

"${git_root_dir}/build/bin/calq-tests"
"${git_root_dir}/build/bin/util-tests"

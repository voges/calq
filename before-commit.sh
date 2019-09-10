#!/usr/bin/env bash

# Apply clang-format in-place on all C++ source code files

set -euo pipefail

git rev-parse --git-dir 1>/dev/null # exit if not inside Git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

{
  echo ""
  echo "==============================================================================="
  echo "${git_root_dir}/ci/cmake-build-release.sh"
  echo "==============================================================================="
  echo ""
  "${git_root_dir}/ci/cmake-build-release.sh";

  echo ""
  echo "==============================================================================="
  echo "${git_root_dir}/ci/cmake-build-debug-with-doc-and-coverage.sh";
  echo "==============================================================================="
  echo ""
  "${git_root_dir}/ci/cmake-build-debug-with-doc-and-coverage.sh";

  echo ""
  echo "==============================================================================="
  echo "${git_root_dir}/ci/run-all-tests.sh";
  echo "==============================================================================="
  echo ""
  "${git_root_dir}/ci/run-all-tests.sh";

  echo ""
  echo "==============================================================================="
  echo "${git_root_dir}/ci/run-all-tests-with-valgrind.sh";
  echo "==============================================================================="
  echo ""
  "${git_root_dir}/ci/run-all-tests-with-valgrind.sh";

  echo ""
  echo "==============================================================================="
  echo "${git_root_dir}/ci/report-code-coverage.sh";
  echo "==============================================================================="
  echo ""
  "${git_root_dir}/ci/report-code-coverage.sh";

  echo ""
  echo "==============================================================================="
  echo "${git_root_dir}/scripts/util/apply-clang-format.sh";
  echo "==============================================================================="
  echo ""
  "${git_root_dir}/scripts/util/apply-clang-format.sh";

  echo ""
  echo "==============================================================================="
  echo "${git_root_dir}/scripts/util/generate-authors-file.sh";
  echo "==============================================================================="
  echo ""
  "${git_root_dir}/scripts/util/generate-authors-file.sh";

} &>"${git_root_dir}/before-commit.sh.log"

echo "Done. Wrote log to: ${git_root_dir}/before-commit.sh.log"

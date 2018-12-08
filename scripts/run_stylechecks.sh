#!/usr/bin/env bash


# -----------------------------------------------------------------------------
# Gather all files to be checked
# -----------------------------------------------------------------------------

# Get the Git root directory
readonly git_root_dir="$(git rev-parse --show-toplevel)"

# Gather all C++ files
cpp_files=()
for node in "$git_root_dir"/source/calq/*; do
    if [ ! -d "$node" ]; then
        cpp_files+=("$node")
    else
        printf "Omitting directory '%s'\\n" "$node"
    fi
done

# Check if all C++ files exist
for cpp_file in "${cpp_files[@]}"; do
    if [ ! -f "$cpp_file" ]; then
        printf "Error: '%s' is not a regular file\\n" "$cpp_file"
        exit -1
    fi
done


# -----------------------------------------------------------------------------
# Check if cppcheck and cpplint.py are available
# -----------------------------------------------------------------------------

# cppcheck
readonly cppcheck="$git_root_dir/external_tools/cppcheck-1.85/build/bin/cppcheck"
if [ ! -f "$cppcheck" ]; then
    printf "Error: '%s' is not a regular file\\n" "$cppcheck"
    exit -1
fi

# cpplint.py
readonly cpplint_py="$git_root_dir/external_tools/cpplint-1.3.0/cpplint.py"
if [ ! -f "$cpplint_py" ]; then
    printf "Error: '%s' is not a regular file\\n" "$cpplint_py"
    exit -1
fi


# -----------------------------------------------------------------------------
# cppcheck
# -----------------------------------------------------------------------------

# Do the work
for cpp_file in "${cpp_files[@]}"; do
    printf "Running cppcheck on: %s\\n" "$cpp_file"
    "$cppcheck" \
        --enable=all \
        --language=c++ \
        --std=posix \
        "$cpp_file"
done


# -----------------------------------------------------------------------------
# cpplint
# -----------------------------------------------------------------------------

# Check for CPPLINT.cfg
cpplint_cfg="$git_root_dir/source/CPPLINT.cfg"
if [ ! -f "$cpplint_cfg" ]; then
    printf "Error: '%s' is not a regular file\\n" "$cpplint_cfg"
    exit -1
fi

# Do the work
for cpp_file in "${cpp_files[@]}"; do
    printf "Running cpplint on: %s\\n" "$cpp_file"
    python "$cpplint_py" "$cpp_file"
done

#!/usr/bin/env bash

# Get the Git root directory
readonly git_root_dir="$(git rev-parse --show-toplevel)"

# We plan to use the 'tmp' directory for our tests - check if it already
# exists; otherwise create it and proceed
readonly tmp_dir="$git_root_dir/tmp"
if [ -d "$tmp_dir" ]; then
    printf "Directory '%s' already exists; aborting calqapp tests\\n" "$tmp_dir"
    exit -1
fi
printf "Creating 'tmp' directory: %s\\n" "$tmp_dir"
mkdir -p "$tmp_dir" || exit -1

# Get the gabacify executable
readonly calqapp="$git_root_dir/build/bin/calqapp"
if [ ! -x "$calqapp" ]; then
    printf "calqapp '%s' is not executable; aborting roundtrip tests\\n" "$calqapp"
    exit -1
fi

# Gather the test files
input_files=()
compressed_files=()
# configuration_files=()
#
input_files+=("$git_root_dir/resources/test_files/fourreads.sam")
compressed_files+=("$git_root_dir/resources/test_files/fourreads.sam.cq")
decompressed_files+=("$git_root_dir/resources/test_files/fourreads.sam.cq.qual")

input_files+=("$git_root_dir/resources/test_files/large_aux.sam")
compressed_files+=("$git_root_dir/resources/test_files/large_aux.sam.cq")
decompressed_files+=("$git_root_dir/resources/test_files/large_aux.sam.cq.qual")


# Do the test roundtrips
for i in "${!input_files[@]}"; do
    printf "Running roundtrip %s\\n" "$i"

    # Gather the i-th files
    input_file=${input_files[$i]}
    compressed_file=${compressed_files[$i]}
    decompressed_file=${decompressed_files[$i]}

    # Encode
    "$calqapp" encode \
        --input_file_path "$input_file" \
        --output_file_path "$compressed_file" 
        || exit -1

    # Decode
    "$calqapp" decode \
        -d \
        --input_file_path "$bytestream_file" \
        --output_file_path "$decompressed_file" \
       || exit -1

    # Check
    head "$decompressed_file" || exit -1

    printf "Success\\n"
    
    rm "$compressed_file" "$decompressed_file"
done

# Delete the 'tmp' directory
rm -rf "$tmp_dir" || exit -1
printf "Removed 'tmp' directory: %s\\n" "$tmp_dir"

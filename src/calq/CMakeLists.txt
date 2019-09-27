project("calq")

set(lib_source_files
        fasta-file-reader.cc
        file-reader.cc
        file-writer.cc
        filter-buffer.cc
        genotyper.cc
        haplotyper.cc
        helpers.cc
        lloyd-max-quantizer.cc
        probability-distribution.cc
        qual-decoder.cc
        qual-encoder.cc
        quantizer.cc
        sam-file-reader.cc
        sam-pileup.cc
        sam-record.cc
        softclip-spreader.cc
        uniform-min-max-quantizer.cc
        uniform-quantizer.cc
)

add_library(calq ${lib_source_files})

set (app_source_files
        calq-codec.cc
        decode.cc
        encode.cc
        program-options.cc
)

add_executable(calq-codec ${app_source_files})

target_link_libraries(calq-codec calq)

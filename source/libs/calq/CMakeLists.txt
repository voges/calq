project("calq")

set(source_files
        errors.h
        exceptions.cc
        filter-buffer.cc
        lloyd-max-quantizer.cc
        probability-distribution.cc
        quantizer.cc
        softclip-spreader.cc
)

add_library(calq ${source_files})

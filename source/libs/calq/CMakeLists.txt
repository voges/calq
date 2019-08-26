project("calq")

set(source_files
        errors.h
        exceptions.cc
        lloyd-max-quantizer.cc
        probability-distribution.cc
        quantizer.cc
)

add_library(calq ${source_files})

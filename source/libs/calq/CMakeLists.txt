project("calq")

set(source_files
        calq-codec.cc
        exceptions.cc
        filter-buffer.cc
        genotyper.cc
        haplotyper.cc
        lloyd-max-quantizer.cc
        log.cc
        probability-distribution.cc
        qual-decoder.cc
        qual-encoder.cc
        quantizer.cc
        sam-pileup.cc
        softclip-spreader.cc
        uniform-min-max-quantizer.cc
        uniform-quantizer.cc
)

add_library(calq ${source_files})

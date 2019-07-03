project("calq")

set(source_files
    calq_coder.cc
    exceptions.cc
    filter_buffer.cc
    genotyper.cc
    haplotyper.cc
    lloyd_max_quantizer.cc
    log.cc
    probability_distribution.cc
    qual_decoder.cc
    qual_encoder.cc
    quantizer.cc
    sam_pileup.cc
    sam_pileup_deque.cc
    softclip_spreader.cc
    uniform_min_max_quantizer.cc
    uniform_quantizer.cc
)

add_library(calq ${source_files})

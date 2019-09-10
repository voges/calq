project("calq")

set(source_files
        calq-codec.cc
        exceptions.cc
        fasta-file-reader.cc
        file.cc
        file-reader.cc
        file-writer.cc
        filter-buffer.cc
        genotyper.cc
        haplotyper.cc
        helpers.cc
        lloyd-max-quantizer.cc
        log.cc
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

add_library(calq ${source_files})

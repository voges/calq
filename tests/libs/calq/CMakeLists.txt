project("calq-tests")

set(source_files
        circular-buffer-test.cc
        fasta-file-reader-test.cc
        file-line-reader-test.cc
        file-reader-test.cc
        file-test.cc
        file-writer-test.cc
        filter-buffer-test.cc
        genotyper-test.cc
        haplotyper-test.cc
        helpers.cc
        helpers-test.cc
        lloyd-max-quantizer-test.cc
        main.cc
        probability-distribution-test.cc
        quantizer-test.cc
        sam-file-reader-test.cc
        sam-pileup-test.cc
        sam-record-test.cc
        softclip-spreader-test.cc
        uniform-min-max-quantizer-test.cc
        uniform-quantizer-test.cc
)

add_executable(calq-tests ${source_files})

target_link_libraries(calq-tests PRIVATE gtest_main)
target_link_libraries(calq-tests PRIVATE calq)

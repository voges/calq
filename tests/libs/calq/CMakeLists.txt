project("calq-tests")

set(source_files
        circular-buffer-test.cc
        filter-buffer-test.cc
        haplotyper-test.cc
        lloyd-max-quantizer-test.cc
        main.cc
        softclip-spreader-test.cc
)

add_executable(calq-tests ${source_files})

target_link_libraries(calq-tests PRIVATE gtest_main)
target_link_libraries(calq-tests PRIVATE calq)

project("calq-tests")

set(source_files
        lloyd-max-quantizer-test.cc
        main.cc
)

add_executable(calq-tests ${source_files})

target_link_libraries(calq-tests PRIVATE gtest_main)
target_link_libraries(calq-tests PRIVATE calq)

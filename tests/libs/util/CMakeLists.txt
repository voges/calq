project("util-tests")

set(source_files
        fasta-file-reader-test.cc
        file-reader-test.cc
        file-test.cc
        file-writer-test.cc
        helpers.cc
        helpers-test.cc
        main.cc
        sam-file-reader-test.cc
        sam-record-test.cc
)

add_executable(util-tests ${source_files})

target_link_libraries(util-tests PRIVATE gtest_main)
target_link_libraries(util-tests PRIVATE util)

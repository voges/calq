project("cip")

set(source_files
        cip.cc
        cq-file.cc
        program-options.cc
        sam-block.cc
        # sam-file-handler.cc
)

add_executable(cip ${source_files})

target_link_libraries(cip calq)
target_link_libraries(cip util)

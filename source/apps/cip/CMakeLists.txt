project("cip")

set(source_files
        cip.cc
        # cq-file.cc
        exceptions.cc
        # fasta-file.cc
        # file.cc
        helpers.cc
        program-options.cc
        # sam-block.cc
        # sam-file-handler.cc
        # sam-file.cc
        # sam-record.cc
)

add_executable(cip ${source_files})

# target_link_libraries(cip calq)

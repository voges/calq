project("cip")

set(source_files
        cip.cc
        program-options.cc
)

add_executable(cip ${source_files})

target_link_libraries(cip calq)

project("cip")

    find_package(Boost 1.54.0 REQUIRED COMPONENTS program_options)

        set(source_files calq.cc cq - file.cc fasta -
            file.cc file.cc helpers.cc program - options.cc sam - block.cc sam -
            file - handler.cc sam - file.cc sam - record.cc)

            add_executable(cip ${source_files})

                target_include_directories(cip PRIVATE ${Boost_INCLUDE_DIR})

                    target_link_libraries(cip calq)
                        target_link_libraries(cip ${Boost_LIBRARIES})

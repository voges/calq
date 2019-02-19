# set project
project("calqapp")

#==============================================================================
# gabac
#==============================================================================

configure_file(CMakeLists_gabac.txt.in gabac-download/CMakeLists.txt)
set(BUILD_GABACIFY OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
execute_process(
    COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result2
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gabac-download
)
if(result)
    message(FATAL_ERROR "CMake step for gabac failed: ${result2}")
endif()
execute_process(
    COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result2
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gabac-download
)
if(result2)
    message(FATAL_ERROR "Build step for gabac failed: ${result2}")
endif()

# Add gabac directly to our build. This defines the gabac targets. 
add_subdirectory(
    ${CMAKE_CURRENT_BINARY_DIR}/gabac-src
    ${CMAKE_CURRENT_BINARY_DIR}/gabac-build
    EXCLUDE_FROM_ALL
)

# Set up include directory
set(gabac_include_dir ${CMAKE_CURRENT_BINARY_DIR}/gabac-src/source)

#==============================================================================
# Boost
#==============================================================================

find_package(Boost 1.54.0 REQUIRED COMPONENTS program_options)

#==============================================================================
# calq exe
#==============================================================================

set(calqapp "calqapp")
set(calqapp_source_dir ${PROJECT_SOURCE_DIR})
set(calqapp_header_dir ${PROJECT_SOURCE_DIR})
get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)
set(calqapp_include_dir ${PARENT_DIR})

# Get all source and header files
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/calq.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/cq_file.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/fasta_file.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/helpers.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/program_options.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/sam_block.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/sam_file.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/sam_record.cc)
set(calqapp_source_files ${calqapp_source_files} ${calqapp_source_dir}/SAMFileHandler.cc)

set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/cq_file.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/error_reporter.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/fasta_file.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/helpers.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/logging.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/program_options.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/sam_block.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/sam_file.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/sam_record.h)
set(calqapp_header_files ${calqapp_header_files} ${calqapp_header_dir}/SAMFileHandler.h)


source_group("Source files" FILES ${calqapp_source_files})
source_group("Header files" FILES ${calqapp_header_files})

add_executable(${calqapp} EXCLUDE_FROM_ALL ${calqapp_source_files} ${calqapp_header_files})

target_include_directories(${calqapp} PRIVATE ${calqapp_include_dir})
target_include_directories(${calqapp} PRIVATE ${calq_include_dir})
target_include_directories(${calqapp} PRIVATE ${Boost_INCLUDE_DIR})
message("Gabac Include dir: ${gabac_include_dir}")
target_include_directories(${calqapp} PRIVATE ${gabac_include_dir})
target_include_directories(${calqapp} PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(${calqapp} ${calq} gabac)
target_link_libraries(${calqapp} ${Boost_LIBRARIES})
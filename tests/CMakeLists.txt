# set project
project("calqtests")

#==============================================================================
# Google Test
#==============================================================================

# Download and unpack googletest at configure time
configure_file(CMakeLists_googletest.txt.in googletest-download/CMakeLists.txt)
execute_process(
    COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
)
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(
    COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
)
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines the gtest and gtest_main
# targets.
add_subdirectory(
    ${CMAKE_CURRENT_BINARY_DIR}/googletest-src
    ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
    EXCLUDE_FROM_ALL
)

# The gtest/gtest_main targets carry header search path dependencies
# automatically when using CMake 2.8.11 or later. Otherwise we have to add
# them here ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()

#==============================================================================
# tests
#==============================================================================

# Set the paths
set(calqtests "calqtests")
set(calq_tests_source_dir ${PROJECT_SOURCE_DIR})
set(calq_tests_header_dir ${PROJECT_SOURCE_DIR})
get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)
set(calq_tests_include_dir ${PARENT_DIR})

# List all source files (alphabetically)
set(calq_tests_source_files ${calq_tests_source_files} ${calq_tests_source_dir}/example_test.cpp)
set(calq_tests_source_files ${calq_tests_source_files} ${calq_tests_source_dir}/main.cpp)

# List all header files (alphabetically)
# set(tests_header_files ${tests_header_files} ${tests_header_dir}/whatever.h)

# Group the source and header files
source_group("Tests Source Files" FILES ${calq_tests_source_files})
# source_group("Tests Header Files" FILES ${tests_header_files})

# Set up the target
add_executable(${calqtests} EXCLUDE_FROM_ALL ${calq_tests_source_files})
target_link_libraries(${calqtests} gtest_main)
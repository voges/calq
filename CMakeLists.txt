#==============================================================================
# Generic
#==============================================================================

# We require at least CMake 3.1
cmake_minimum_required(VERSION 3.1)
message(STATUS "CMake generator: ${CMAKE_GENERATOR}")

# C and C++ code shall be compliant to C(++)11
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

# This is calq
project("calq")


#==============================================================================
# Compiler
#==============================================================================

message(STATUS "CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
else()
    message(FATAL_ERROR "C++ compiler currently not supported by CMake build")
endif()


#==============================================================================
# Build type
#==============================================================================

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Build type not specified")
    set(CMAKE_BUILD_TYPE Debug)
    message(STATUS "Using fall-back build type '${CMAKE_BUILD_TYPE}'")
endif()

# Log the compiler flags for each build type
message(STATUS "C++ compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Additional C++ compiler flags for build type 'Debug': ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "Additional C++ compiler flags for build type 'Release': ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Additional C++ compiler flags for build type 'RelWithDebInfo': ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "Additional C++ compiler flags for build type 'MinSizeRel': ${CMAKE_CXX_FLAGS_MINSIZEREL}")
message(STATUS "Selected build type: ${CMAKE_BUILD_TYPE}")


#==============================================================================
# Google Test
#==============================================================================

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
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
# calq
#==============================================================================

# Set the paths
set(calq "calq")
set(calq_source_dir ${CMAKE_SOURCE_DIR}/source/calq)
set(calq_header_dir ${CMAKE_SOURCE_DIR}/source/calq)
set(calq_include_dir ${CMAKE_SOURCE_DIR}/source)

# Get all source and header files
file(GLOB_RECURSE calq_c_source_files ${calq_source_dir}/*.c)
file(GLOB_RECURSE calq_cc_source_files ${calq_source_dir}/*.cc)
file(GLOB_RECURSE calq_header_files ${calq_header_dir}/*.h)

# Group the source and header files
source_group("CALQ Source Files" FILES ${calq_c_source_files} ${calq_cc_source_files})
source_group("CALQ Header Files" FILES ${calq_header_files})

# Set up the target
add_executable(${calq} ${calq_c_source_files} ${calq_cc_source_files})
target_include_directories(${calq} PRIVATE ${calq_include_dir})
target_include_directories(${calq} PRIVATE ${calq_include_dir}/calq)


#==============================================================================
# tests
#==============================================================================

# Set the paths
set(tests "tests")
set(tests_source_dir ${CMAKE_SOURCE_DIR}/tests)
set(tests_header_dir ${CMAKE_SOURCE_DIR}/tests)
set(tests_include_dir ${CMAKE_SOURCE_DIR})

# List all source files (alphabetically)
set(tests_source_files ${tests_source_files} ${tests_source_dir}/example_test.cpp)
set(tests_source_files ${tests_source_files} ${tests_source_dir}/main.cpp)

# List all header files (alphabetically)
# set(tests_header_files ${tests_header_files} ${tests_header_dir}/whatever.h)

# Group the source and header files
source_group("Tests Source Files" FILES ${tests_source_files})
# source_group("Tests Header Files" FILES ${tests_header_files})

# Set up the target
add_executable(${tests} EXCLUDE_FROM_ALL ${tests_source_files})
target_link_libraries(${tests} gtest_main)

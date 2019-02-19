#==============================================================================
# Generic
#==============================================================================

# we require at least CMake 3.1
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# C and C++ code shall be compliant to C(++)11
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

# Initializing our build options
option(CALQ_BUILD_LIB "Builds the calq library" ON)
option(CALQ_BUILD_APP "Builds the calq library and application" ON)
option(CALQ_BUILD_TESTS "Builds the googletests for calq" ON)

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

# Set output dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set build dirs
if(CALQ_BUILD_LIB OR CALQ_BUILD_APP)
	add_subdirectory(source/calq)
endif()

if(CALQ_BUILD_APP)
	add_subdirectory(source/calqapp)
endif()

if(CALQ_BUILD_TESTS AND CALQ_BUILD_APP)
	add_subdirectory(tests/)
endif()

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

# Set up project-specific variables
set(PROJECT_NAME "calq")
project(${PROJECT_NAME})
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_ROOT_DIR ${CMAKE_SOURCE_DIR})
set(PROJECT_BUILD_DIR ${CMAKE_BINARY_DIR})
set(PROJECT_SOURCE_DIR ${PROJECT_ROOT_DIR}/src/)
set(PROJECT_INCLUDE_DIR ${PROJECT_ROOT_DIR}/src/)

# Get all source and header files (C and C++)
file(GLOB_RECURSE C_SOURCE_FILES ${PROJECT_SOURCE_DIR}/*.c)
file(GLOB_RECURSE CC_SOURCE_FILES ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB_RECURSE CPP_SOURCE_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE H_HEADER_FILES ${PROJECT_INCLUDE_DIR}/*.h)
file(GLOB_RECURSE HH_HEADER_FILES ${PROJECT_INCLUDE_DIR}/*.hh)
file(GLOB_RECURSE HPP_HEADER_FILES ${PROJECT_INCLUDE_DIR}/*.hpp)
set(PROJECT_SOURCE_FILES ${C_SOURCE_FILES} ${CC_SOURCE_FILES} ${CPP_SOURCE_FILES})
set(PROJECT_HEADER_FILES ${H_HEADER_FILES} ${HH_HEADER_FILES} ${HPP_HEADER_FILES})

# Group the source and header files
source_group("Source files" FILES ${PROJECT_SOURCE_FILES})
source_group("Header files" FILES ${PROJECT_HEADER_FILES})

# Set compiler flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Detected GNU Compiler Collection")

    # C compiler flags
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -02")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-equal")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wpointer-arith")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-align")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wstrict-prototypes")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wwrite-strings")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Waggregate-return")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-qual")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wswitch-default")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wswitch-enum")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wunreachable-code")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wconversion")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wundef")

    # C++ compiler flags
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(STATUS "Detected (Apple)Clang")

    # C compiler flags
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -02")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")

    # C++ compiler flags
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(STATUS "Detected Microsoft Visual Studio Compiler")
else ()
    message(FATAL_ERROR "Compiler not supported: ${CMAKE_CXX_COMPILER_ID}")
endif ()

# Generate version.h
file(WRITE ${PROJECT_BUILD_DIR}/version.cmake "file(WRITE version.h \"#ifndef CALQ_VERSION_H_\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define CALQ_VERSION_H_\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define CALQ_VERSION \\\"${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}\\\"\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define CALQ_VERSION_MAJOR ${PROJECT_VERSION_MAJOR}\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define CALQ_VERSION_MINOR ${PROJECT_VERSION_MINOR}\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define CALQ_VERSION_PATCH ${PROJECT_VERSION_PATCH}\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#endif // CALQ_VERSION_H_\\n\\n\")\n")
add_custom_target(version
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BUILD_DIR}/version.cmake
    ADD_DEPENDENCIES ${PROJECT_BUILD_DIR}/version.cmake
)

# Write build time in UTC ISO 8601 to timestamp.h
message(STATUS "Creating target for ${PROJECT_BUILD_DIR}/timestamp.h")
file(WRITE ${PROJECT_BUILD_DIR}/timestamp.cmake "STRING(TIMESTAMP TIMESTAMP_UTC UTC)\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(WRITE timestamp.h \"#ifndef CALQ_TIMESTAMP_H_\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(APPEND timestamp.h \"#define CALQ_TIMESTAMP_H_\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(APPEND timestamp.h \"#define CALQ_TIMESTAMP_UTC \\\"\${TIMESTAMP_UTC}\\\"\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(APPEND timestamp.h \"#endif // CALQ_TIMESTAMP_H_\\n\\n\")\n")
add_custom_target(timestamp
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BUILD_DIR}/timestamp.cmake
    ADD_DEPENDENCIES ${PROJECT_BUILD_DIR}/timestamp.cmake
)

# Generate gitrevision.h if Git is available and the .git directory is present
find_program(GIT_SCM git DOC "Git version control")
mark_as_advanced(GIT_SCM)
find_file(GIT_DIR NAMES .git PATHS ${PROJECT_ROOT_DIR} NO_DEFAULT_PATH)
if (GIT_SCM AND EXISTS ${GIT_DIR})

    # Create git.h that depends on the Git HEAD log
    message(STATUS "Git is available and .git/ directory is present")
    message(STATUS "Creating target for ${PROJECT_BUILD_DIR}/git.h")

    # Get the current working branch
    execute_process(COMMAND ${GIT_SCM} rev-parse --abbrev-ref HEAD
                    WORKING_DIRECTORY ${PROJECT_ROOT_DIR}
                    OUTPUT_VARIABLE GIT_BRANCH
                    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Get the latest abbreviated commit hash of the working branch
    execute_process(COMMAND ${GIT_SCM} log -1 --format=%h #"--pretty=format:%H"
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_COMMIT_HASH_SHORT
                    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Get the latest long commit hash of the working branch
    execute_process(COMMAND ${GIT_SCM} log -1 --format=%H #"--pretty=format:%H"
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_COMMIT_HASH_LONG
                    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    file(WRITE ${PROJECT_BUILD_DIR}/git.cmake "file(WRITE git.h \"#ifndef CALQ_GIT_H_\\n\")\n")
    file(APPEND ${PROJECT_BUILD_DIR}/git.cmake "file(APPEND git.h \"#define CALQ_GIT_H_\\n\\n\")\n")
    file(APPEND ${PROJECT_BUILD_DIR}/git.cmake "file(APPEND git.h \"#define CALQ_GIT_BRANCH \\\"${GIT_BRANCH}\\\"\\n\")\n")
    file(APPEND ${PROJECT_BUILD_DIR}/git.cmake "file(APPEND git.h \"#define CALQ_GIT_COMMIT_HASH_SHORT \\\"${GIT_COMMIT_HASH_SHORT}\\\"\\n\")\n")
    file(APPEND ${PROJECT_BUILD_DIR}/git.cmake "file(APPEND git.h \"#define CALQ_GIT_COMMIT_HASH_LONG \\\"${GIT_COMMIT_HASH_LONG}\\\"\\n\\n\")\n")
    file(APPEND ${PROJECT_BUILD_DIR}/git.cmake "file(APPEND git.h \"#endif // CALQ_GIT_H_\\n\\n\")\n")
    add_custom_target(git
        COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BUILD_DIR}/git.cmake
        ADD_DEPENDENCIES ${GIT_DIR}/logs/HEAD
    )
else ()
    # No version control
    message(STATUS "Git is not available or .git/ directory not found")
    message(FATAL_ERROR "Build stopped")
endif ()

# Doxygen documentation
find_package(Doxygen)
set(DOXYFILE_IN ${PROJECT_ROOT_DIR}/Doxyfile.in)
if (DOXYGEN_FOUND AND EXISTS ${DOXYFILE_IN})
    message(STATUS "Doxygen is available")
    message(STATUS "Found Doxyfile.in: ${DOXYFILE_IN}")
    message(STATUS "Issue 'make doc' to build documentation")
    set(DOXYFILE ${PROJECT_BUILD_DIR}/Doxyfile)
    set(DOXYGEN_DIR ${PROJECT_BUILD_DIR}/doc/)
    configure_file(${DOXYFILE_IN} ${DOXYFILE} @ONLY)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
        WORKING_DIRECTORY ${PROJECT_BUILD_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
else ()
    message(WARNING "Doxygen or Doxyfile.in not found; target 'doc' not available")
endif ()

# Includes, targets, and dependencies
include_directories(${PROJECT_BUILD_DIR})
include_directories(${PROJECT_INCLUDE_DIR})
add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILES} ${PROJECT_HEADER_FILES})
add_dependencies(${PROJECT_NAME} timestamp)
add_dependencies(${PROJECT_NAME} git)
#add_dependencies(${PROJECT_NAME} doc)
add_dependencies(${PROJECT_NAME} version)
#target_link_libraries(${PROJECT_NAME} z)


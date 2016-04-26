cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

# set up project-specific variables
set(PROJECT_NAME "calq")
project(${PROJECT_NAME})
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_ROOT_DIR ${CMAKE_SOURCE_DIR})
set(PROJECT_BUILD_DIR ${CMAKE_BINARY_DIR})
set(PROJECT_SOURCE_DIR ${PROJECT_ROOT_DIR}/src)

# get all source and header files (C and C++)
file(GLOB_RECURSE C_SOURCE_FILES ${PROJECT_SOURCE_DIR}/*.c)
file(GLOB_RECURSE CC_SOURCE_FILES ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB_RECURSE CPP_SOURCE_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE H_HEADER_FILES ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB_RECURSE HH_HEADER_FILES ${PROJECT_SOURCE_DIR}/*.hh)
file(GLOB_RECURSE HPP_HEADER_FILES ${PROJECT_SOURCE_DIR}/*.hpp)
set(PROJECT_SOURCE_FILES ${C_SOURCE_FILES} ${CC_SOURCE_FILES} ${CPP_SOURCE_FILES})
set(PROJECT_HEADER_FILES ${H_HEADER_FILES} ${HH_HEADER_FILES} ${HPP_HEADER_FILES})

# C compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -02")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-equal")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wpointer-arith")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-align")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wstrict-prototypes")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wwrite-strings")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Waggregate-return")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-qual")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wswitch-default")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wswitch-enum")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wunreachable-code")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wconversion")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wundef")

# CXX compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")

# generate version.h file
message(STATUS "Creating file ${PROJECT_BUILD_DIR}/version.h")
file(WRITE ${PROJECT_BUILD_DIR}/version.cmake "file(WRITE version.h \"#ifndef CMAKE_VERSION_H\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define CMAKE_VERSION_H\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define _CMAKE_VERSION_MAJOR_ ${PROJECT_VERSION_MAJOR}\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define _CMAKE_VERSION_MINOR_ ${PROJECT_VERSION_MINOR}\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#define _CMAKE_VERSION_PATCH_ ${PROJECT_VERSION_PATCH}\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/version.cmake "file(APPEND version.h \"#endif // CMAKE_VERSION_H\\n\\n\")\n")
add_custom_target(version
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BUILD_DIR}/version.cmake
    ADD_DEPENDENCIES ${PROJECT_BUILD_DIR}/version.cmake
)

# write build time in UTC ISO 8601 to timestamp.h
message(STATUS "Creating file ${PROJECT_BUILD_DIR}/timestamp.h")
file(WRITE ${PROJECT_BUILD_DIR}/timestamp.cmake "STRING(TIMESTAMP UTCTIME UTC)\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(WRITE timestamp.h \"#ifndef CMAKE_TIMESTAMP_H\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(APPEND timestamp.h \"#define CMAKE_TIMESTAMP_H\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(APPEND timestamp.h \"#define _CMAKE_UTCTIME_ \\\"\${UTCTIME}\\\"\\n\\n\")\n")
file(APPEND ${PROJECT_BUILD_DIR}/timestamp.cmake "file(APPEND timestamp.h \"#endif // CMAKE_TIMESTAMP_H\\n\\n\")\n")
add_custom_target(timestamp
    COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BUILD_DIR}/timestamp.cmake
    ADD_DEPENDENCIES ${PROJECT_BUILD_DIR}/timestamp.cmake
)

# generate gitrevision.h if Git is available and the .git directory is present
find_program(GIT_SCM git DOC "Git version control")
mark_as_advanced(GIT_SCM)
find_file(GIT_DIR NAMES .git PATHS ${PROJECT_ROOT_DIR} NO_DEFAULT_PATH)
if (GIT_SCM AND GIT_DIR)
    # create gitrevision.h that depends on the Git HEAD log
    message(STATUS "Git is available and found .git/ directory")
    message(STATUS "Creating file ${PROJECT_BUILD_DIR}/gitrevision.h")

    add_custom_command(
        OUTPUT ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "#ifndef CMAKE_GITREVISION_H" > ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "#define CMAKE_GITREVISION_H" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo_append "#define _CMAKE_GITREVISION_LONG_ \"" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${GIT_SCM} log -1 "--pretty=format:%H" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "\"" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo_append "#define _CMAKE_GITREVISION_SHORT_ \"" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${GIT_SCM} log -1 "--pretty=format:%h" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "\"" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "#endif // CMAKE_GITREVISION_H" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${PROJECT_BUILD_DIR}/gitrevision.h
        DEPENDS ${GIT_DIR}/logs/HEAD
        VERBATIM
    )
else ()
    # no version control
    message(STATUS "Will not remake ${PROJECT_BUILD_DIR}/gitrevision.h")
endif ()
add_custom_target(gitrevision DEPENDS ${PROJECT_BUILD_DIR}/gitrevision.h)

# Doxygen documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
    message(STATUS "Doxygen is available; issue 'make doc' to build documentation")
    set(DOXYGEN_DOXYFILE_IN ${PROJECT_ROOT_DIR}/Doxyfile.in)
    set(DOXYGEN_DOXYFILE_OUT ${PROJECT_BUILD_DIR}/Doxyfile)
    set(DOXYGEN_OUTPUT_DIR ${PROJECT_BUILD_DIR}/doc/)
    configure_file(${DOXYGEN_DOXYFILE_IN} ${DOXYGEN_DOXYFILE_OUT} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_DOXYFILE_OUT}
        WORKING_DIRECTORY ${PROJECT_BUILD_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
else ()
    message(WARNING "Doxygen not found")
endif ()

# includes, targets, and dependencies
include_directories(${PROJECT_BUILD_DIR})
include_directories(${PROJECT_SOURCE_DIR})
add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_FILES})
add_dependencies(${PROJECT_NAME} timestamp)
add_dependencies(${PROJECT_NAME} gitrevision)
#add_dependencies(${PROJECT_NAME} doc)
add_dependencies(${PROJECT_NAME} version)


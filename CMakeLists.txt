cmake_minimum_required(VERSION 2.6)

# set up project-specific variables
set(PROJ_NAME "calq")
project(${PROJ_NAME})
set(PROJ_SOURCE_DIR ${CMAKE_SOURCE_DIR})
set(PROJ_BINARY_DIR ${CMAKE_BINARY_DIR})
set(PROJ_INCLUDES ${PROJ_SOURCE_DIR}/src)
set(PROJ_INCLUDES ${PROJ_INCLUDES} ${PROJ_SOURCE_DIR}/src/qualcodec)

# get all source and header files
file(GLOB_RECURSE PROJ_SOURCES src/*.c)
file(GLOB_RECURSE PROJ_HEADERS src/*.h)

# compiler flags
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -Wall)
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

# comment to disable debug output
add_definitions("-DDBG")

# write build time in UTC ISO 8601 to timestamp.h
message(STATUS "Creating file ${PROJ_BINARY_DIR}/timestamp.h")
file(WRITE ${PROJ_BINARY_DIR}/timestamp.cmake "STRING(TIMESTAMP UTCTIME UTC)\n")
file(APPEND ${PROJ_BINARY_DIR}/timestamp.cmake "FILE(WRITE timestamp.h \"#ifndef CMAKE_TIMESTAMP_H\\n\")\n")
file(APPEND ${PROJ_BINARY_DIR}/timestamp.cmake "FILE(APPEND timestamp.h \"#define CMAKE_TIMESTAMP_H\\n\\n\")\n")
file(APPEND ${PROJ_BINARY_DIR}/timestamp.cmake "FILE(APPEND timestamp.h \"#define _CMAKE_UTCTIME_ \\\"\${UTCTIME}\\\"\\n\\n\")\n")
file(APPEND ${PROJ_BINARY_DIR}/timestamp.cmake "FILE(APPEND timestamp.h \"#endif // CMAKE_TIMESTAMP_H\\n\\n\")\n")
add_custom_target (
    timestamp
    COMMAND ${CMAKE_COMMAND} -P ${PROJ_BINARY_DIR}/timestamp.cmake
    ADD_DEPENDENCIES ${PROJ_BINARY_DIR}/timestamp.cmake
)

# generate gitrevision.h if Git is available and the .git directory is present
find_program(GIT_SCM git DOC "Git version control")
mark_as_advanced(GIT_SCM)
find_file(GITDIR NAMES .git PATHS ${CMAKE_SOURCE_DIR} NO_DEFAULT_PATH)
if (GIT_SCM AND GITDIR)
    # create gitrevision.h that depends on the Git HEAD log
    message(STATUS "Git is available")
    message(STATUS "Found directory .git/")
    message(STATUS "Creating file ${PROJ_BINARY_DIR}/gitrevision.h")

    add_custom_command(
        OUTPUT ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "#ifndef CMAKE_GITREVISION_H" > ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "#define CMAKE_GITREVISION_H" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo_append "#define _CMAKE_GITREVISION_LONG_ \"" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${GIT_SCM} log -1 "--pretty=format:%H" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "\"" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo_append "#define _CMAKE_GITREVISION_SHORT_ \"" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${GIT_SCM} log -1 "--pretty=format:%h" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "\"" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "#endif // CMAKE_GITREVISION_H" >> ${PROJ_BINARY_DIR}/gitrevision.h
        COMMAND ${CMAKE_COMMAND} -E echo "" >> ${PROJ_BINARY_DIR}/gitrevision.h
        DEPENDS ${GITDIR}/logs/HEAD
        VERBATIM
    )
else()
    # no version control
    message(STATUS "Will not remake ${PROJ_BINARY_DIR}/gitrevision.h")
endif()
add_custom_target(gitrevision DEPENDS ${PROJ_BINARY_DIR}/gitrevision.h)

# includes, targets, and dependencies
include_directories(${PROJ_BINARY_DIR})
include_directories(${PROJ_INCLUDES})
add_executable(${PROJ_NAME} ${PROJ_SOURCES})
add_dependencies(${PROJ_NAME} timestamp)
add_dependencies(${PROJ_NAME} gitrevision)


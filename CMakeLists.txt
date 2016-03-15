cmake_minimum_required(VERSION 2.6)

set(PROJ_NAME "calq")
set(PROJ_SOURCE_PATH ${CMAKE_SOURCE_DIR})
set(PROJ_BIN_PATH ${CMAKE_BINARY_DIR})
set(PROJ_INCLUDES ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/qualcodec)

file(GLOB_RECURSE PROJ_SOURCES src/*.c)
file(GLOB_RECURSE PROJ_HEADERS src/*.h)

project(${PROJ_NAME})

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-equal")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wpointer-arith")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-align")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wstrict-prototypes")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wwrite-strings")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Waggregate-return")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-qual")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wswitch-default")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wswitch-enum")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wunreachable-code")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wconversion")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wundef")

add_definitions("-DDBG") # comment to disable debug output

include_directories(${PROJ_INCLUDES})
add_executable(${PROJ_NAME} ${PROJ_SOURCES})


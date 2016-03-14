cmake_minimum_required(VERSION 2.6)
project(calq)

set (CMAKE_CXX_FLAGS "-D DBG")
include_directories(. qualcodec)
file(GLOB SOURCES "*.c" "qualcodec/*.c")
add_executable(calq ${SOURCES})
install(TARGETS calq RUNTIME DESTINATION bin)

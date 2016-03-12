cmake_minimum_required(VERSION 2.6)
project(calq)

add_executable(calq main.cpp)

install(TARGETS calq RUNTIME DESTINATION bin)

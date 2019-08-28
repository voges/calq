#==============================================================================
# Generic
#==============================================================================

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project("calq")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "Using fall-back build type '${CMAKE_BUILD_TYPE}'")
endif()

include(AddGoogleTest)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)


#==============================================================================
# Compiler
#==============================================================================

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")


#==============================================================================
# Executable/library paths
#==============================================================================

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


#==============================================================================
# Include paths
#==============================================================================

include_directories("${CMAKE_SOURCE_DIR}/source/apps/")
include_directories("${CMAKE_SOURCE_DIR}/source/libs/")
include_directories("${CMAKE_SOURCE_DIR}/source/libs/third-party")


#==============================================================================
# Subdirectories
#==============================================================================

add_subdirectory(source)
add_subdirectory(tests)


#==============================================================================
# Doxygen documentation
#==============================================================================

#find_package(Doxygen)
#if (DOXYGEN_FOUND)
#    message(STATUS "Doxygen is available")
#    message(STATUS "Issue 'make doc' to build documentation")
#
#    if (NOT DOXYGEN_DOT_FOUND)
#        message(WARNING "!!!! Graphviz's `dot' tool not found on this system. "
#                "Doxygen documentation will not include graphs.")
#    endif ()
#
#    # For forward compatibility with newer versions of cmake
#    if(NOT DEFINED DOXYGEN_HAVE_DOT)
#        set(DOXYGEN_HAVE_DOT ${DOXYGEN_DOT_FOUND})
#    endif ()
#
#    set(DOXYGEN_STRIP_FROM_PATH      ${CMAKE_SOURCE_DIR})
#    set(DOXYGEN_ENUM_VALUES_PER_LINE 1)
#    # set input and output files
#    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
#    set(DOXYGEN_OUT ${MPEGG_DOC_DIR}/Doxyfile)
#    # request to configure the file
#    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#
#    add_custom_target(doc
#            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#            WORKING_DIRECTORY ${MPEGG_DOC_DIR}
#            COMMENT "Generating API documentation with Doxygen"
#            VERBATIM)
#else ()
#    message(WARNING "Doxygen not found; target 'doc' not available")
#endif ()

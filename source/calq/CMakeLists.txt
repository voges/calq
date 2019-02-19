#==============================================================================
# calq
#==============================================================================

# set project
project("calq")

# Set the paths
set(calq "calq")
set(calq_source_dir ${PROJECT_SOURCE_DIR})
set(calq_header_dir ${PROJECT_SOURCE_DIR})
get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)
set(calq_include_dir ${PARENT_DIR})

# Get all source and header files
set(calq_source_files ${calq_source_files} ${calq_source_dir}/calq_coder.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/exceptions.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/filter_buffer.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/genotyper.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/haplotyper.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/lloyd_max_quantizer.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/log.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/probability_distribution.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/qual_decoder.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/qual_encoder.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/sam_pileup.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/sam_pileup_deque.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/softclip_spreader.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/uniform_min_max_quantizer.cc)
set(calq_source_files ${calq_source_files} ${calq_source_dir}/uniform_quantizer.cc)


set(calq_header_files ${calq_header_files} ${calq_header_dir}/calq.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/calq_coder.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/circular_buffer.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/error_exception_reporter.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/exceptions.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/filter_buffer.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/genotyper.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/haplotyper.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/lloyd_max_quantizer.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/log.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/probability_distribution.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/qual_decoder.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/qual_encoder.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/sam_pileup.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/sam_pileup_deque.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/softclip_spreader.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/uniform_min_max_quantizer.h)
set(calq_header_files ${calq_header_files} ${calq_header_dir}/uniform_quantizer.h)

# Group the source and header files
source_group("CALQ Source Files" FILES ${calq_source_files})
source_group("CALQ Header Files" FILES ${calq_header_files})

set(CALQ_BUILD_SHARED_LIB OFF CACHE BOOL "Build shared CALQ library.")
if (${CALQ_BUILD_SHARED_LIB})
    message(STATUS "Configured build of CALQ as *shared* library")
    add_library(${calq} SHARED ${calq_source_files} ${calq_header_files})
    set_target_properties(${calq} PROPERTIES LINKER_LANGUAGE CXX)
    set_target_properties(${calq} PROPERTIES POSITION_INDEPENDENT_CODE 1)
else ()
    message(STATUS "Configured build of calq as *static* library")
    add_library(${calq} STATIC ${calq_source_files} ${calq_header_files})
    set_target_properties(${calq} PROPERTIES LINKER_LANGUAGE CXX)
endif ()

target_include_directories(${calq} PRIVATE ${calq_include_dir})
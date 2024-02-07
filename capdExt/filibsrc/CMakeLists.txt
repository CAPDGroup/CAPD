cmake_minimum_required(VERSION 3.13)
project(filib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

###############################################################################################################################
# BUILD LIBRARY
###############################################################################################################################

add_library(${PROJECT_NAME})
target_compile_options(${PROJECT_NAME} PUBLIC -O2 -frounding-math)
target_compile_options(${PROJECT_NAME} PRIVATE -DFILIB_EXTENDED)

if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
    message(STATUS "Processor arch: 64-bit x86 (SSE)")
    target_compile_options(${PROJECT_NAME} PUBLIC -DFILIB_HAVE_SSE)
else()
    message(FATAL_ERROR "Unknown or unsupported processor architecture.")
endif()

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ieee/primitive.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/interval/stdfun/filib_consts.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_double_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_float_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_base_double_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_base_float_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_sse_const.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_x87_const.cpp)

###############################################################################################################################
# BUILD TEST CASES
###############################################################################################################################

if ( ${CAPD_BUILD_TESTS} )

    enable_testing()
    add_subdirectory(test)

endif()

###############################################################################################################################
# INSTALLATION
###############################################################################################################################

set(HEADER_FILES)
file(GLOB_RECURSE HEADER_FILES_HPP ${CMAKE_CURRENT_SOURCE_DIR}/ieee/*.hpp)
list(APPEND HEADER_FILES ${HEADER_FILES_HPP})

file(GLOB_RECURSE HEADER_FILES_HPP ${CMAKE_CURRENT_SOURCE_DIR}/interval/*.hpp)
list(APPEND HEADER_FILES ${HEADER_FILES_HPP})

file(GLOB_RECURSE HEADER_FILES_HPP ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/*.hpp)
list(APPEND HEADER_FILES ${HEADER_FILES_HPP})

file(GLOB_RECURSE HEADER_FILES_HPP ${CMAKE_CURRENT_SOURCE_DIR}/rounding_control/*.hpp)
list(APPEND HEADER_FILES ${HEADER_FILES_HPP})

file(GLOB_RECURSE HEADER_FILES_ICC ${CMAKE_CURRENT_SOURCE_DIR}/*.icc)
list(APPEND HEADER_FILES ${HEADER_FILES_ICC})

foreach(HEADER_FILE ${HEADER_FILES})
    file(RELATIVE_PATH TARGET_DIR ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE})
    string(REGEX MATCH "(.*\/)" DIR ${TARGET_DIR})

    install(FILES ${HEADER_FILE} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${DIR}")
endforeach()

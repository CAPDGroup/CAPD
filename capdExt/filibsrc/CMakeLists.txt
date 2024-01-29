cmake_minimum_required(VERSION 3.13)
project(filib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

###############################################################################################################################
# BUILD LIBRARY
###############################################################################################################################

add_library(${PROJECT_NAME})
target_compile_options(${PROJECT_NAME} PUBLIC -O0 -fno-fast-math -frounding-math)
target_compile_options(${PROJECT_NAME} PRIVATE -DFILIB_EXTENDED)

if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
    message(STATUS "Processor arch: 64-bit x86 (SSE)")
    target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_SSE)
else()
    message(FATAL_ERROR "Unknown or unsupported processor architecture.")
endif()

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

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

enable_testing()
add_subdirectory(test)







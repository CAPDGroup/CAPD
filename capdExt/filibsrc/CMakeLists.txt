#cmake_minimum_required(VERSION 3.10)
#project(libprim)
#
#add_library(libprim)

target_include_directories(${PROJECT_NAME} PUBLIC .)

target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ieee/primitive.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/interval/stdfun/filib_consts.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_double_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_float_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_base_double_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_base_float_setup.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_sse_const.cpp )
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fp_traits/fp_traits_x87_const.cpp)
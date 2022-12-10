cmake_minimum_required(VERSION 3.10)
project(capd_experimental)

message("Project " ${CMAKE_PROJECT_NAME} " configuration started")

add_library(${CMAKE_PROJECT_NAME})

enable_testing()

add_subdirectory(capdExt)
add_subdirectory(capdAux)
add_subdirectory(capdAlg)
add_subdirectory(capdDynSys4)

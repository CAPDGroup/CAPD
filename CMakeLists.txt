cmake_minimum_required(VERSION 3.10)
project(capd)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/capdMake/cmake/")
include(CAPDVersion.txt)
include(CAPD_options)
include(CAPD_functions)

add_library(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

target_compile_options(${PROJECT_NAME} PUBLIC -O2 -fno-fast-math -frounding-math)

enable_testing()


add_subdirectory(capdExt)
add_subdirectory(capdAux)
add_subdirectory(capdAlg)
add_subdirectory(capdDynSys4)

if (${CAPD_ENABLE_MULTIPRECISION})

    target_link_libraries(${PROJECT_NAME} PUBLIC -lmpfr)
    target_link_libraries(${PROJECT_NAME} PUBLIC -lgmp)

    target_compile_options(${PROJECT_NAME} PUBLIC -D__HAVE_MPFR__)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION "lib")

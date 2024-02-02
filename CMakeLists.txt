cmake_minimum_required(VERSION 3.13)
project(capd)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/capdMake/cmake/")
include(CAPDVersion.txt)
include(CAPD_options)
include(CAPD_functions)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

add_library(${PROJECT_NAME})
add_dependencies(${PROJECT_NAME} filib)
target_compile_options(${PROJECT_NAME} PUBLIC -D__USE_FILIB__)

target_compile_options(${PROJECT_NAME} PUBLIC -O2 -frounding-math)
target_link_libraries(${PROJECT_NAME} PUBLIC filib)

enable_testing()

if (${CAPD_ENABLE_MULTIPRECISION})

    target_link_libraries(${PROJECT_NAME} PUBLIC -lmpfr)
    target_link_libraries(${PROJECT_NAME} PUBLIC -lgmp)

    target_compile_options(${PROJECT_NAME} PUBLIC -D__HAVE_MPFR__)
endif()

add_subdirectory(capdExt)
add_subdirectory(capdAux)
add_subdirectory(capdAlg)
add_subdirectory(capdDynSys)
add_subdirectory(capdMake/libcapd)
add_subdirectory(capdMake/bin)

install(TARGETS
        ${PROJECT_NAME} filib
    EXPORT CapdConfig
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    INCLUDES DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
)

export(TARGETS
        ${PROJECT_NAME} filib
    NAMESPACE Capd::
    FILE "${CMAKE_CURRENT_BINARY_DIR}/CapdConfig.cmake"
)
install(EXPORT
        CapdConfig
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Capd"
    NAMESPACE Capd::
)

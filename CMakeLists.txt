cmake_minimum_required(VERSION 3.13)
project(capd)

message(STATUS "CMake version: ${CMAKE_VERSION}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/capdMake/cmake/")
include(CAPDVersion.txt)
include(CAPD_options)
include(CAPD_functions)

if (${CAPD_INTERVAL_TYPE} STREQUAL "TIGHT")
    set(CMAKE_CXX_STANDARD 20) # required by the TIGHT library
else()
    set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

add_library(${PROJECT_NAME})
target_compile_options(${PROJECT_NAME} PUBLIC -O2 -frounding-math)

if (${CAPD_INTERVAL_TYPE} STREQUAL "NATIVE")

    target_compile_options(${PROJECT_NAME} PUBLIC -D__USE_NATIVE__)

elseif(${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")

    add_dependencies(${PROJECT_NAME} filib)
    target_compile_options(${PROJECT_NAME} PUBLIC -D__USE_FILIB__)
    target_link_libraries(${PROJECT_NAME} PUBLIC filib)


elseif(${CAPD_INTERVAL_TYPE} STREQUAL "TIGHT")

    # TIGHT: correctly rounded intervals; it fetches its own dependencies (NFG, CORE-MATH) with CPM
    if(NOT EXISTS "${CAPD_TIGHT_DIR}/CMakeLists.txt")
        find_package(Git QUIET)
        if(NOT GIT_FOUND)
            message(FATAL_ERROR "CAPD_INTERVAL_TYPE=TIGHT: git is needed to download TIGHT (or set CAPD_TIGHT_DIR to its sources)")
        endif()
        # 1. TIGHT is a submodule of CAPD (capdExt/tight): fetch it together with its submodules
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git" AND "${CAPD_TIGHT_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/capdExt/tight")
            message(STATUS "Downloading TIGHT (git submodule capdExt/tight)")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- capdExt/tight
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} RESULT_VARIABLE _tight_res)
        endif()
        # 2. CAPD is not a git checkout (e.g. a source archive): clone the pinned revision
        if(NOT EXISTS "${CAPD_TIGHT_DIR}/CMakeLists.txt")
            message(STATUS "Downloading TIGHT (git clone ${CAPD_TIGHT_GIT_URL})")
            file(REMOVE_RECURSE "${CAPD_TIGHT_DIR}")
            execute_process(COMMAND ${GIT_EXECUTABLE} clone --recurse-submodules ${CAPD_TIGHT_GIT_URL} ${CAPD_TIGHT_DIR} RESULT_VARIABLE _tight_res)
            if(_tight_res EQUAL 0)
                execute_process(COMMAND ${GIT_EXECUTABLE} checkout ${CAPD_TIGHT_GIT_TAG} WORKING_DIRECTORY ${CAPD_TIGHT_DIR} RESULT_VARIABLE _tight_res)
                execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive WORKING_DIRECTORY ${CAPD_TIGHT_DIR})
            endif()
        endif()
        if(NOT EXISTS "${CAPD_TIGHT_DIR}/CMakeLists.txt")
            message(FATAL_ERROR "Could not download TIGHT from ${CAPD_TIGHT_GIT_URL}; clone it manually into ${CAPD_TIGHT_DIR} or set CAPD_TIGHT_DIR")
        endif()
    endif()
    if(NOT DEFINED ENV{CPM_SOURCE_CACHE})
        set(ENV{CPM_SOURCE_CACHE} "${CAPD_TIGHT_DIR}/.cpm")
    endif()
    set(TIGHT_CORE_MATH ON CACHE BOOL "" FORCE)
    add_subdirectory(${CAPD_TIGHT_DIR} ${CMAKE_BINARY_DIR}/tight EXCLUDE_FROM_ALL)
    add_dependencies(${PROJECT_NAME} tight)
    # tight is built with -mfma; do not let the compiler fuse a*b+c in CAPD code (results of plain double computations would differ from other builds)
    target_compile_options(${PROJECT_NAME} PUBLIC -D__USE_TIGHT__ -DTIGHT_CORE_MATH -ffp-contract=off)
    target_link_libraries(${PROJECT_NAME} PUBLIC tight::tight)

elseif(${CAPD_INTERVAL_TYPE} STREQUAL "CXSC")

    add_dependencies(${PROJECT_NAME} cxsc)
    target_compile_options(${PROJECT_NAME} PUBLIC -D__USE_CXSC__)
    target_link_libraries(${PROJECT_NAME} PUBLIC cxsc)

else()

    message(FATAL_ERROR "Unsupported CAPD_INTERVAL_TYPE")

endif()

enable_testing()

if (${CAPD_ENABLE_MULTIPRECISION})

    find_library(MPFR_LIB_FOUND mpfr)
    find_library(GMP_LIB_FOUND gmp)

    if (MPFR_LIB_FOUND)
        message(STATUS "Library `mpfr` found.")
    else()
        message(FATAL_ERROR "Library `mpfr` not found. Please install it or set CAPD_ENABLE_MULTIPRECISION to OFF.")
    endif()

    if (GMP_LIB_FOUND)
        message(STATUS "Library `gmp` found.")
    else()
        message(FATAL_ERROR "Library `gmp` not found. Please install it or set CAPD_ENABLE_MULTIPRECISION to OFF.")
    endif()

    target_link_libraries(${PROJECT_NAME} PUBLIC -lmpfr)
    target_link_libraries(${PROJECT_NAME} PUBLIC -lgmp)

    target_compile_options(${PROJECT_NAME} PUBLIC -D__HAVE_MPFR__)
endif()

add_subdirectory(capdExt)
add_subdirectory(capdAux)
add_subdirectory(capdAlg)
add_subdirectory(capdDynSys)
add_subdirectory(capdMake/libcapd)

if (${CAPD_INTERVAL_TYPE} STREQUAL "NATIVE")

    install(TARGETS
            ${PROJECT_NAME}
        EXPORT capd-config
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        INCLUDES DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
    )

    export(TARGETS
            ${PROJECT_NAME}
        NAMESPACE capd::
        FILE "${CMAKE_CURRENT_BINARY_DIR}/capd-config.cmake"
    )


elseif(${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")

    install(TARGETS
            ${PROJECT_NAME} filib
        EXPORT capd-config
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        INCLUDES DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
    )

    export(TARGETS
            ${PROJECT_NAME} filib
        NAMESPACE capd::
        FILE "${CMAKE_CURRENT_BINARY_DIR}/capd-config.cmake"
    )

elseif (${CAPD_INTERVAL_TYPE} STREQUAL "TIGHT")

    # TIGHT and its dependencies are built in-tree (static); no CMake package is exported.
    # Use bin/capd-config (pkg-config) to get flags for user programs.
    message(STATUS "CAPD with TIGHT: capd-config.cmake is not exported, use bin/capd-config")

elseif (${CAPD_INTERVAL_TYPE} STREQUAL "CXSC")

    install(TARGETS
            ${PROJECT_NAME} cxsc
        EXPORT capd-config
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        INCLUDES DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
    )

    export(TARGETS
            ${PROJECT_NAME} cxsc
        NAMESPACE capd::
        FILE "${CMAKE_CURRENT_BINARY_DIR}/capd-config.cmake"
    )

else()

    message(FATAL_ERROR "Unsupported CAPD_INTERVAL_TYPE")

endif()

if (NOT "${CAPD_INTERVAL_TYPE}" STREQUAL "TIGHT")
install(EXPORT
        capd-config
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/capd"
        NAMESPACE capd::
    )
endif()

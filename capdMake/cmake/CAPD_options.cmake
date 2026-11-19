option(CAPD_BUILD_ALL "CAPD test and examples " OFF)
# message(STATUS "Build all CAPD      : ${CAPD_BUILD_ALL}")
option(CAPD_BUILD_TESTS "Activate/deactivate CAPD tests compilation" OFF)
option(CAPD_BUILD_EXAMPLES "Activate/deactivate CAPD examples compilation" OFF)

if( ${CAPD_BUILD_ALL} )
  set(CAPD_BUILD_TESTS ON)
  set(CAPD_BUILD_EXAMPLES ON)
endif()

option(CAPD_ENABLE_MULTIPRECISION "Activate/decactivate multiprecision support in CAPD" OFF)

set(CAPD_INTERVAL_TYPE "FILIB" CACHE STRING "Select interval type (NATIVE, FILIB, CXSC, TIGHT)")

set_property(CACHE CAPD_INTERVAL_TYPE PROPERTY STRINGS "NATIVE" "FILIB" "CXSC" "TIGHT")

# TIGHT intervals (correctly rounded, https://gitlab.com/fsichetti/tight): path to its source directory
set(CAPD_TIGHT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../tight" CACHE PATH "Source directory of the TIGHT library (used when CAPD_INTERVAL_TYPE=TIGHT)")

message(STATUS "Build CAPD tests    : ${CAPD_BUILD_TESTS}")
message(STATUS "Build CAPD examples : ${CAPD_BUILD_EXAMPLES}")
message(STATUS "CAPD multiprecision : ${CAPD_ENABLE_MULTIPRECISION}")
message(STATUS "CAPD interval type  : ${CAPD_INTERVAL_TYPE}")

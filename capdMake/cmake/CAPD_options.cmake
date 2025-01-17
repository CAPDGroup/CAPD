option(CAPD_BUILD_ALL "CAPD test and examples " OFF)
# message(STATUS "Build all CAPD      : ${CAPD_BUILD_ALL}")
option(CAPD_BUILD_TESTS "Activate/deactivate CAPD tests compilation" OFF)
option(CAPD_BUILD_EXAMPLES "Activate/deactivate CAPD examples compilation" OFF)

if( ${CAPD_BUILD_ALL} )
  set(CAPD_BUILD_TESTS ON)
  set(CAPD_BUILD_EXAMPLES ON)
endif()

option(CAPD_ENABLE_MULTIPRECISION "Activate/decactivate multiprecision support in CAPD" OFF)

set(CAPD_INTERVAL_TYPE "FILIB" CACHE STRING "Select interval type (NATIVE, FILIB, CXSC)")

set_property(CACHE CAPD_INTERVAL_TYPE PROPERTY STRINGS "NATIVE" "FILIB" "CXSC")

message(STATUS "Build CAPD tests    : ${CAPD_BUILD_TESTS}")
message(STATUS "Build CAPD examples : ${CAPD_BUILD_EXAMPLES}")
message(STATUS "CAPD multiprecision : ${CAPD_ENABLE_MULTIPRECISION}")
message(STATUS "CAPD interval type  : ${CAPD_INTERVAL_TYPE}")

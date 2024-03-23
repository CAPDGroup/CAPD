option(CAPD_BUILD_ALL "CAPD test and examples " OFF)
# message(STATUS "Build all CAPD      : ${CAPD_BUILD_ALL}")
option(CAPD_BUILD_TESTS "Activate/deactivate CAPD tests compilation" OFF)
option(CAPD_BUILD_EXAMPLES "Activate/deactivate CAPD examples compilation" OFF)

if( ${CAPD_BUILD_ALL} )
  set(CAPD_BUILD_TESTS ON)
  set(CAPD_BUILD_EXAMPLES ON)
endif()

option(CAPD_ENABLE_MULTIPRECISION "Activate/decactivate multiprecision support in CAPD" OFF)

message(STATUS "Build CAPD tests    : ${CAPD_BUILD_TESTS}")
message(STATUS "Build CAPD examples : ${CAPD_BUILD_EXAMPLES}")
message(STATUS "CAPD multiprecision : ${CAPD_ENABLE_MULTIPRECISION}")

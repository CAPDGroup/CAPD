# option(BUILD_TEST_EXECUTABLES "Build test executables" OFF)
# option(BUILD_EXAMPLE_EXECUTABLES "Build example executables" OFF)
# option(BUILD_KRAK "Build krak component" OFF)
# option(ENABLE_MULTIPRECISION "Enable multiprecision (requires mpfr and gmp)" OFF)


#option(WITH_CAPD_TESTS "Activate/deactivate tests compilation" OFF)
#option(WITH_CAPD_EXAMPLES "Activate/deactivate tests compilation" OFF)

#if( ${WITH_CAPD_TESTS} )
#  set(BUILD_TEST_EXECUTABLES ON)
#else ()
#  set(BUILD_TEST_EXECUTABLES OFF)
#endif ()
option(CAPD_BUILD_ALL "CAPD test and examples " OFF)
message(STATUS "Build all CAPD      : ${CAPD_BUILD_ALL}")
option(BUILD_TEST_EXECUTABLES "Activate/deactivate CAPD tests compilation" OFF)
option(BUILD_EXAMPLES_EXECUTABLES "Activate/deactivate CAPD examples compilation" OFF)

#cmake_dependent_option(BUILD_TEST_EXECUTABLES "Activate/deactivate CAPD tests compilation" OFF "NOT CAPD_BUILD_ALL" ON )
if( ${CAPD_BUILD_ALL} )
  message(STATUS "TAK")
  set(BUILD_TEST_EXECUTABLES ON)
  set(BUILD_EXAMPLES_EXECUTABLES ON)
endif()

option(ENABLE_MULTIPRECISION "Activate/decactivate multiprecision support in CAPD" ON)

message(STATUS "Build CAPD tests    : ${BUILD_TEST_EXECUTABLES}")
message(STATUS "Build CAPD examples : ${BUILD_EXAMPLES_EXECUTABLES}")
message(STATUS "CAPD multiprecision : ${ENABLE_MULTIPRECISION}")
#  if ( ${BUILD_EXAMPLE_EXECUTABLES} ))

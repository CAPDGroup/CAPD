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

# TIGHT intervals (correctly rounded, https://gitlab.com/minimize-solve/tight), used when CAPD_INTERVAL_TYPE=TIGHT.
# TIGHT is a git submodule in capdExt/tight; it is downloaded automatically at the configuration step if it is missing
# (together with its dependencies NFG, SIMDe and CORE-MATH, which TIGHT itself fetches into capdExt/tight/.cpm).
set(CAPD_TIGHT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/capdExt/tight" CACHE PATH "Source directory of the TIGHT library (used when CAPD_INTERVAL_TYPE=TIGHT)")
set(CAPD_TIGHT_GIT_URL "https://gitlab.com/minimize-solve/tight.git" CACHE STRING "Git repository of TIGHT (used when the submodule is not available)")
set(CAPD_TIGHT_GIT_TAG "486875d57dee62beec8aeb25e5aa7e21e668bce3" CACHE STRING "Revision of TIGHT used when it is cloned without the submodule")

message(STATUS "Build CAPD tests    : ${CAPD_BUILD_TESTS}")
message(STATUS "Build CAPD examples : ${CAPD_BUILD_EXAMPLES}")
message(STATUS "CAPD multiprecision : ${CAPD_ENABLE_MULTIPRECISION}")
message(STATUS "CAPD interval type  : ${CAPD_INTERVAL_TYPE}")

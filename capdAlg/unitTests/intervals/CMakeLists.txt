# Comparison of the interval arithmetics (FILIB, TIGHT, CXSC, NATIVE) - correctness (versus MPFR) and speed.
# The backend is selected by CAPD_INTERVAL_TYPE, so every backend needs its own build directory;
# compare_backends.py puts the results of several builds side by side.
if (${CAPD_ENABLE_MULTIPRECISION})

    add_executable(IntervalCompare intervalCompare.cpp)
    target_link_libraries(IntervalCompare capd)

    # short version (a few seconds); the full version: IntervalCompare --full
    add_test(NAME IntervalCompareQuick COMMAND IntervalCompare --quick --csv ${CMAKE_CURRENT_BINARY_DIR}/intervals_${CAPD_INTERVAL_TYPE}.csv)

    # quality of enclosures for point arguments and arguments of diameter 1e-13, 1e-10, 1e-8
    # (full version: IntervalCompare --full --width); --strict: fails if a result does not contain the exact hull
    add_test(NAME IntervalCompareWidthQuick COMMAND IntervalCompare --quick --width --strict --csv ${CMAKE_CURRENT_BINARY_DIR}/intervals_width_${CAPD_INTERVAL_TYPE}.csv)

endif()

# The old gtest based tests (they are not maintained and need the gtest library)
option(CAPD_BUILD_LEGACY_INTERVAL_GTESTS "Build the old gtest based tests of intervals" OFF)
if (CAPD_BUILD_LEGACY_INTERVAL_GTESTS)

    set(TEST_PROJECT_NAME IntervalTest)
    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_INSTALL_PREFIX}/lib)
    target_link_libraries(${TEST_PROJECT_NAME} capd)

    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DIntervalTest.cpp)
    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DoubleIntervalTest.cpp)
    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DoubleRoundingTest.cpp)
    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IntervalTest.cpp)

    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endif()

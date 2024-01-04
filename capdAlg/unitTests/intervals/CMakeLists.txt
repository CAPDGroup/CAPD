set(TEST_PROJECT_NAME IntervalTest)
add_executable(${TEST_PROJECT_NAME})
target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(${TEST_PROJECT_NAME} capd)

target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DIntervalTest.cpp)
target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DoubleIntervalTest.cpp)
target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DoubleRoundingTest.cpp)
target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IntervalTest.cpp)

add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

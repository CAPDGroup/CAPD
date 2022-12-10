set(TEST_PROJECT_NAME C2CurveMaskTest)
add_executable(${TEST_PROJECT_NAME})
target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)
find_package (Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)
target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/C2CurveMaskTest.cpp)

add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

set(TEST_PROJECT_NAME C1CurveMaskTest)
add_executable(${TEST_PROJECT_NAME})
target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)
find_package (Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)
target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/C1CurveMaskTest.cpp)

add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

set(TEST_PROJECT_NAME CnCurveMaskTest)
add_executable(${TEST_PROJECT_NAME})
target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)
find_package (Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)
target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/CnCurveMaskTest.cpp)

add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

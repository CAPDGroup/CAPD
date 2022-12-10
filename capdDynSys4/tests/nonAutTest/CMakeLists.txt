set(TEST_PROJECT_NAME nonAutTest)
add_executable(${TEST_PROJECT_NAME})
target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)
target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/nonAutTest.cpp)

add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

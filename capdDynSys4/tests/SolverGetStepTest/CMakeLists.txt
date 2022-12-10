set (TEST_PROJECT_NAME SolverGetStepTest)
message("${TEST_PROJECT_NAME}")

add_executable(${TEST_PROJECT_NAME})

target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)

find_package (Boost REQUIRED COMPONENTS unit_test_framework)
# target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${Boost_LIBRARY_DIRS})
target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)

target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IFixStep.cpp)
# target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/SolverGetStepTest.cpp)
# target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IVarStep.cpp)
# target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DFixStep.cpp)
# target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/DVarStep.cpp)

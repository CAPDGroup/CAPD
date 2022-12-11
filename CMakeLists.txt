cmake_minimum_required(VERSION 3.10)
project(capd_experimental)

message("Project " ${CMAKE_PROJECT_NAME} " configuration started")

add_library(${CMAKE_PROJECT_NAME})

enable_testing()

function(build_example_executable_monofile EXAMPLE_EXECUTABLE_NAME)

    add_executable(${EXAMPLE_EXECUTABLE_NAME})
    target_link_directories(${EXAMPLE_EXECUTABLE_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
    target_link_libraries(${EXAMPLE_EXECUTABLE_NAME} capd_experimental)
    target_sources(${EXAMPLE_EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_EXECUTABLE_NAME}.cpp)

endfunction()

function(build_test_executable_monofile TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
    target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)
    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PROJECT_NAME}.cpp)

    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

function(build_boost_test_executable_monofile TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
    target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)
    find_package (Boost REQUIRED COMPONENTS unit_test_framework)
    target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)
    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PROJECT_NAME}.cpp)
    
    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

function(build_test_executable TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
    target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)

    foreach(SOURCE_FILE IN LISTS ARGN)
        target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_FILE})
    endforeach()

    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

function(build_boost_test_executable TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
    target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)
    find_package (Boost REQUIRED COMPONENTS unit_test_framework)
    target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)

    foreach(SOURCE_FILE IN LISTS ARGN)
        target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_FILE})
    endforeach()

    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

set(BUILD_TEST_EXECUTABLES true)
set(BUILD_EXAMPLE_EXECUTABLES true)

add_subdirectory(capdExt)
add_subdirectory(capdAux)
add_subdirectory(capdAlg)
add_subdirectory(capdDynSys4)

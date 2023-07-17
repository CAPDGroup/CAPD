cmake_minimum_required(VERSION 3.10)
project(capd2)

set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

add_library(${PROJECT_NAME})

target_compile_options(${PROJECT_NAME} PUBLIC -O1 -fno-fast-math -frounding-math)

enable_testing()

function(build_example_executable_monofile EXAMPLE_EXECUTABLE_NAME)

    add_executable(${EXAMPLE_EXECUTABLE_NAME})
    target_link_directories(${EXAMPLE_EXECUTABLE_NAME} PRIVATE ${CMAKE_INSTALL_PREFIX}/lib)
    target_link_libraries(${EXAMPLE_EXECUTABLE_NAME} capd2)
    target_sources(${EXAMPLE_EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_EXECUTABLE_NAME}.cpp)

endfunction()

function(build_test_executable_monofile TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_INSTALL_PREFIX}/lib)
    target_link_libraries (${TEST_PROJECT_NAME} capd2)
    target_link_libraries (${TEST_PROJECT_NAME} -lstdc++)
    target_link_libraries (${TEST_PROJECT_NAME} -lm)
    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PROJECT_NAME}.cpp)

    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

function(build_boost_test_executable_monofile TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_INSTALL_PREFIX}/lib)
    target_link_libraries(${TEST_PROJECT_NAME} capd2)
    find_package (Boost REQUIRED COMPONENTS unit_test_framework)
    target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)
    target_link_libraries (${TEST_PROJECT_NAME} -lstdc++)
    target_link_libraries (${TEST_PROJECT_NAME} -lm)
    target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PROJECT_NAME}.cpp)
    
    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

function(build_test_executable TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_INSTALL_PREFIX}/lib)
    target_link_libraries (${TEST_PROJECT_NAME} capd2)
    target_link_libraries (${TEST_PROJECT_NAME} -lstdc++)
    target_link_libraries (${TEST_PROJECT_NAME} -lm)

    foreach(SOURCE_FILE IN LISTS ARGN)
        target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_FILE})
    endforeach()

    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

function(build_boost_test_executable TEST_PROJECT_NAME)

    add_executable(${TEST_PROJECT_NAME})
    target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_INSTALL_PREFIX}/lib)
    target_link_libraries(${TEST_PROJECT_NAME} capd2)
    find_package (Boost REQUIRED COMPONENTS unit_test_framework)
    target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)
    target_link_libraries (${TEST_PROJECT_NAME} -lstdc++)
    target_link_libraries (${TEST_PROJECT_NAME} -lm)

    foreach(SOURCE_FILE IN LISTS ARGN)
        target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_FILE})
    endforeach()

    add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

function(install_headers)

    foreach(HEADER_FILE IN LISTS ARGN)
        string(REGEX MATCH "(.*\/)" DIR ${HEADER_FILE})
        install(FILES ${HEADER_FILE} DESTINATION "include/${DIR}")
    endforeach()

endfunction()

add_subdirectory(capdExt)
add_subdirectory(capdAux)
add_subdirectory(capdAlg)
add_subdirectory(capdDynSys4)

if (${ENABLE_MULTIPRECISION})
    target_link_libraries(${PROJECT_NAME} PUBLIC -lmpfr)
    target_link_libraries(${PROJECT_NAME} PUBLIC -lgmp)

    target_compile_options(${PROJECT_NAME} PUBLIC -D__HAVE_MPFR__)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION "lib")

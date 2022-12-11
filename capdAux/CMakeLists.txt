target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC include)

add_subdirectory(src/capd/auxil)
add_subdirectory(src/capd/krak)

if ( ${BUILD_TEST_EXECUTABLES} )

    add_subdirectory(tests/auxilTests)

endif()


if ( ${BUILD_EXAMPLE_EXECUTABLES} )

    add_subdirectory(examples/auxildemo)
    # add_subdirectory(examples/configFileReader)
    # add_subdirectory(examples/krakdemo)
    # add_subdirectory(examples/logger)

endif()

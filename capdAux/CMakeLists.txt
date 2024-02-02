target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

add_subdirectory(include)
add_subdirectory(src/capd/auxil)

if ( ${CAPD_BUILD_TESTS} )

    add_subdirectory(tests/auxilTests)

endif()


if ( ${CAPD_BUILD_EXAMPLES} )

    add_subdirectory(examples/auxildemo)
    # add_subdirectory(examples/configFileReader)
    # add_subdirectory(examples/krakdemo)
    # add_subdirectory(examples/logger)

endif()

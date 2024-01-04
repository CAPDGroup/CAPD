target_include_directories(${PROJECT_NAME} PUBLIC include)
add_subdirectory(include)

add_subdirectory(src/capd/basicalg)
add_subdirectory(src/capd/intervals)
add_subdirectory(src/capd/matrixAlgorithms)
add_subdirectory(src/capd/rings)
add_subdirectory(src/capd/rounding)
add_subdirectory(src/capd/vectalg)

if ( ${ENABLE_MULTIPRECISION} )

    add_subdirectory(src/mpcapd/intervals)
    add_subdirectory(src/mpcapd/matrixAlgorithms)
    add_subdirectory(src/mpcapd/vectalg)

endif()

if ( ${BUILD_TEST_EXECUTABLES} )

    add_subdirectory(tests)
    # add_subdirectory(unitTests/intervals)            # disabled due to build failures

endif()

if ( ${BUILD_EXAMPLE_EXECUTABLES} )

    add_subdirectory(examples/alglibdemo)
    add_subdirectory(examples/intervals)
    add_subdirectory(examples/matrixAlgorithms)
    add_subdirectory(examples/TexWriter)
    add_subdirectory(examples/vectalg)

    if ( ${ENABLE_MULTIPRECISION} )
    
        # add_subdirectory(examples/BigIntegers)              # disabled due to build failures

    endif()

endif()

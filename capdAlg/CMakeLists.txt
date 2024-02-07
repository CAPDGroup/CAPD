target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

add_subdirectory(include)

add_subdirectory(src/capd/basicalg)
add_subdirectory(src/capd/intervals)
add_subdirectory(src/capd/matrixAlgorithms)
add_subdirectory(src/capd/rings)
add_subdirectory(src/capd/rounding)
add_subdirectory(src/capd/vectalg)

if ( ${CAPD_ENABLE_MULTIPRECISION} )

    add_subdirectory(src/mpcapd/intervals)
    add_subdirectory(src/mpcapd/matrixAlgorithms)
    add_subdirectory(src/mpcapd/vectalg)

endif()

if ( ${CAPD_BUILD_TESTS} )

    add_subdirectory(tests)
    # add_subdirectory(unitTests/intervals)            # disabled due to build failures

endif()

if ( ${CAPD_BUILD_EXAMPLES} )

    add_subdirectory(examples/alglibdemo)
    add_subdirectory(examples/intervals)
    add_subdirectory(examples/matrixAlgorithms)
    add_subdirectory(examples/TexWriter)
    add_subdirectory(examples/vectalg)

    if ( ${CAPD_ENABLE_MULTIPRECISION} )
    
        # add_subdirectory(examples/BigIntegers)              # disabled due to build failures

    endif()

endif()

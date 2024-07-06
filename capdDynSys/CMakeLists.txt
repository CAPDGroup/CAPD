target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)


add_subdirectory(include)

add_subdirectory(src/capd/autodiff)
add_subdirectory(src/capd/covrel)
add_subdirectory(src/capd/diffAlgebra)
add_subdirectory(src/capd/diffIncl)
add_subdirectory(src/capd/dynset)
add_subdirectory(src/capd/dynsys)
add_subdirectory(src/capd/geomset)
add_subdirectory(src/capd/map)
add_subdirectory(src/capd/pdes)
add_subdirectory(src/capd/poincare)

if ( ${CAPD_ENABLE_MULTIPRECISION} )

    add_subdirectory(src/mpcapd/autodiff)
    add_subdirectory(src/mpcapd/diffAlgebra)
    add_subdirectory(src/mpcapd/diffIncl)
    add_subdirectory(src/mpcapd/dynset)
    add_subdirectory(src/mpcapd/dynsys)
    add_subdirectory(src/mpcapd/geomset)
    add_subdirectory(src/mpcapd/map)
    add_subdirectory(src/mpcapd/poincare)

endif()

if ( ${CAPD_BUILD_TESTS} )

    add_subdirectory(tests)

endif()


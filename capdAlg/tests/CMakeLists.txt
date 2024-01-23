add_subdirectory(basicalg)
add_subdirectory(complex)
add_subdirectory(intervalTest)
add_subdirectory(intvtst)
add_subdirectory(matrixAlgorithms)
add_subdirectory(vecttst)

if ( ${CAPD_ENABLE_MULTIPRECISION} )

    add_subdirectory(mptests)

endif()

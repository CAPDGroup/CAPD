target_include_directories(${PROJECT_NAME} PUBLIC include)
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

# add_subdirectory(src/mpcapd/autodiff)
# add_subdirectory(src/mpcapd/diffAlgebra)
# add_subdirectory(src/mpcapd/diffIncl)
# add_subdirectory(src/mpcapd/dynset)
# add_subdirectory(src/mpcapd/dynsys)
# add_subdirectory(src/mpcapd/geomset)
# add_subdirectory(src/mpcapd/map)
# add_subdirectory(src/mpcapd/poincare)

if ( ${BUILD_TEST_EXECUTABLES} )

    add_subdirectory(tests)

endif()

if ( ${BUILD_EXAMPLE_EXECUTABLES} )

    add_subdirectory(examples/cndemo)
    add_subdirectory(examples/coveringRelation)
    add_subdirectory(examples/diffIncldemo)
    add_subdirectory(examples/encloseTrajectory)
    # add_subdirectory(examples/graphicalSet)          # disabled due to krak dependency
    add_subdirectory(examples/integrate)
    # add_subdirectory(examples/invariantSet)          # disabled due to krak dependency
    add_subdirectory(examples/LorenzPeriodicOrbit)
    add_subdirectory(examples/maps)
    # add_subdirectory(examples/multiPrecision)        # disabled due to mprecision dependency
    add_subdirectory(examples/newton)
    add_subdirectory(examples/odes)
    add_subdirectory(examples/odesrig)
    add_subdirectory(examples/poincare)
    # add_subdirectory(examples/pointst)               # disabled due to krak dependency
    add_subdirectory(examples/RosslerChaoticDynamics)
    add_subdirectory(examples/RosslerPeriodicOrbit)
    # add_subdirectory(examples/taylor)                # disabled due to krak dependency
    
endif()

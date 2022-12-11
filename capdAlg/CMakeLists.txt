target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC include)

add_subdirectory(src/capd/basicalg)
add_subdirectory(src/capd/intervals)
add_subdirectory(src/capd/matrixAlgorithms)
add_subdirectory(src/capd/rings)
add_subdirectory(src/capd/rounding)
add_subdirectory(src/capd/vectalg)

add_subdirectory(src/mpcapd/intervals)
add_subdirectory(src/mpcapd/matrixAlgorithms)
add_subdirectory(src/mpcapd/vectalg)

add_subdirectory(tests)

# -----------------------------------------------------------------------------------------------------------------------------
# Build example executables
# -----------------------------------------------------------------------------------------------------------------------------
if (true)

    add_subdirectory(examples/alglibdemo)
    # add_subdirectory(examples/BigIntegers)           # disabled due to mprecision dependency
    add_subdirectory(examples/intervals)
    add_subdirectory(examples/matrixAlgorithms)
    add_subdirectory(examples/TexWriter)
    add_subdirectory(examples/vectalg)

endif()



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

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

add_subdirectory(include)

add_subdirectory(src/capd/alglib)
add_subdirectory(filibsrc)


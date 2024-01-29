target_include_directories(${PROJECT_NAME} PUBLIC include)
add_subdirectory(include)

add_subdirectory(src/capd/alglib)
add_subdirectory(filibsrc)


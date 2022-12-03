message("Entering: " ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC include)

add_subdirectory(src/capd/auxil)
add_subdirectory(src/capd/krak)

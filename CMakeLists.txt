cmake_minimum_required(VERSION 3.10)
project(capd_experimental)

message("Project " ${CMAKE_PROJECT_NAME} " configuration started")

add_library(${CMAKE_PROJECT_NAME})

enable_testing()

function(build_example_executable SOURCE_FILE)

    set(EXAMPLE_EXECUTABLE_NAME ${SOURCE_FILE})
    add_executable(${EXAMPLE_EXECUTABLE_NAME})
    target_link_directories(${EXAMPLE_EXECUTABLE_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
    target_link_libraries(${EXAMPLE_EXECUTABLE_NAME} capd_experimental)
    target_sources(${EXAMPLE_EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_FILE}.cpp)

endfunction()

add_subdirectory(capdExt)
add_subdirectory(capdAux)
add_subdirectory(capdAlg)
add_subdirectory(capdDynSys4)

set(EXAMPLE_EXECUTABLE_NAME cndemo1)
add_executable(${EXAMPLE_EXECUTABLE_NAME})
target_link_directories(${EXAMPLE_EXECUTABLE_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${EXAMPLE_EXECUTABLE_NAME} capd_experimental)
target_sources(${EXAMPLE_EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cndemo1.cpp)


set(EXAMPLE_EXECUTABLE_NAME cndemo2)
add_executable(${EXAMPLE_EXECUTABLE_NAME})
target_link_directories(${EXAMPLE_EXECUTABLE_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${EXAMPLE_EXECUTABLE_NAME} capd_experimental)
target_sources(${EXAMPLE_EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cndemo2.cpp)


set(EXAMPLE_EXECUTABLE_NAME cndemo3)
add_executable(${EXAMPLE_EXECUTABLE_NAME})
target_link_directories(${EXAMPLE_EXECUTABLE_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
target_link_libraries(${EXAMPLE_EXECUTABLE_NAME} capd_experimental)
target_sources(${EXAMPLE_EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cndemo3.cpp)

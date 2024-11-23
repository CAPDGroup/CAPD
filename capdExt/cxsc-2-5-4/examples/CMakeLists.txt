
include_directories(${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/CToolbox/Modules
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/fi_lib
        ${CMAKE_SOURCE_DIR}/src/rts)

link_directories(
  ${CMAKE_SOURCE_DIR/src}
)

add_executable(allzeros allzeros.cpp)
add_executable(example example.cpp)
add_executable(inewton inewton.cpp)
add_executable(io io.cpp)
add_executable(lexample lexample.cpp)
add_executable(linewton linewton.cpp)
add_executable(rungekutta rungekutta.cpp)
add_executable(trace trace.cpp)

target_link_libraries(allzeros cxsc)
target_link_libraries(example cxsc)
target_link_libraries(inewton cxsc)
target_link_libraries(io cxsc)
target_link_libraries(lexample cxsc)
target_link_libraries(linewton cxsc)
target_link_libraries(rungekutta cxsc)
target_link_libraries(trace cxsc)

file(GLOB example_src "*.cpp")
install(FILES ${example_src} DESTINATION examples)
install(TARGETS allzeros example inewton io lexample linewton rungekutta trace RUNTIME DESTINATION examples)



include_directories(${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/CToolbox/Modules
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/fi_lib
        ${CMAKE_SOURCE_DIR}/src/rts)

link_directories(
  ${CMAKE_SOURCE_DIR/src}
)

add_executable(hess_test hess_test.cpp )
target_link_libraries(hess_test cxsc)

add_executable(cpz_test cpz_test.cpp )
target_link_libraries(cpz_test cxsc)

add_executable(gop1_test gop1_test.cpp )
target_link_libraries(gop1_test cxsc)

add_executable(gop_test gop_test.cpp )
target_link_libraries(gop_test cxsc)

add_executable(jac_test jac_test.cpp )
target_link_libraries(jac_test cxsc)

add_executable(lop_test lop_test.cpp )
target_link_libraries(lop_test cxsc)

add_executable(lss_test lss_test.cpp )
target_link_libraries(lss_test cxsc)

add_executable(nlfz_test nlfz_test.cpp )
target_link_libraries(nlfz_test cxsc)

add_executable(nlss_test nlss_test.cpp )
target_link_libraries(nlss_test cxsc)

add_executable(rpe_test rpe_test.cpp )
target_link_libraries(rpe_test cxsc)

add_executable(xev_test1 xev_test1.cpp)
target_link_libraries(xev_test1 cxsc)

add_executable(xev_test2 xev_test2.cpp)
target_link_libraries(xev_test2 cxsc)

add_executable(xev_test3 xev_test3.cpp)
target_link_libraries(xev_test3 cxsc)

#added 2013-01-22
file(GLOB ctoolbox_programms_src "*.cpp")
#install(FILES ${ctoolbox_programms_src} DESTINATION testcases)
#install(TARGETS hess_test cpz_test gop1_test gop_test jac_test lop_test lss_test nlfz_test nlss_test rpe_test xev_test1 xev_test2 xev_test3 RUNTIME DESTINATION testcases)

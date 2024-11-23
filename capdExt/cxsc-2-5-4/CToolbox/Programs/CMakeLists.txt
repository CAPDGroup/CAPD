
include_directories(${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/CToolbox/Modules
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/fi_lib
        ${CMAKE_SOURCE_DIR}/src/rts)

link_directories(
   ${CMAKE_INSTALL_PREFIX}/lib
)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(hess_ex hess_ex.cpp)
target_link_libraries(hess_ex cxsc)

add_executable(cpz_ex cpz_ex.cpp)
target_link_libraries(cpz_ex cxsc)

add_executable(ddf_ex ddf_ex.cpp)
target_link_libraries(ddf_ex cxsc)

add_executable(fastlss_ex fastlss_ex.cpp)
target_link_libraries(fastlss_ex cxsc)

add_executable(gop1_ex gop1_ex.cpp)
target_link_libraries(gop1_ex cxsc)

add_executable(gop_ex gop_ex.cpp)
target_link_libraries(gop_ex cxsc)

add_executable(jac_ex jac_ex.cpp)
target_link_libraries(jac_ex cxsc)

add_executable(lop_ex lop_ex.cpp)
target_link_libraries(lop_ex cxsc)

add_executable(lss_ex lss_ex.cpp)
target_link_libraries(lss_ex cxsc)

add_executable(nlfz_ex nlfz_ex.cpp)
target_link_libraries(nlfz_ex cxsc)

add_executable(nlss_ex nlss_ex.cpp)
target_link_libraries(nlss_ex cxsc)

add_executable(rpe_ex rpe_ex.cpp)
target_link_libraries(rpe_ex cxsc)

add_executable(xev_ex1 xev_ex1.cpp)
target_link_libraries(xev_ex1 cxsc)

add_executable(xev_ex2 xev_ex2.cpp)
target_link_libraries(xev_ex2 cxsc)

add_executable(xev_ex3 xev_ex3.cpp)
target_link_libraries(xev_ex3 cxsc)

#added 2012-12-11
file(GLOB ctoolbox_programms_src "*.cpp")
install(FILES ${ctoolbox_programms_src} DESTINATION examples)
install(TARGETS hess_ex cpz_ex ddf_ex fastlss_ex gop1_ex gop_ex jac_ex lop_ex lss_ex nlfz_ex nlss_ex rpe_ex xev_ex1 xev_ex2 xev_ex3 RUNTIME DESTINATION examples)

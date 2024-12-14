
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/CToolbox/Modules
		${CMAKE_CURRENT_SOURCE_DIR}/CToolbox/test
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/fi_lib
        ${CMAKE_CURRENT_SOURCE_DIR}/src/rts)

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

add_executable(tester tester.cpp)
target_link_libraries(tester cxsc)

add_executable(sumup sumup.cpp)
target_link_libraries(sumup cxsc)


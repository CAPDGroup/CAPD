target_include_directories(${PROJECT_NAME} PUBLIC include)
add_subdirectory(include)

add_subdirectory(src/capd/alglib)

include(ExternalProject)
ExternalProject_Add(
    filib
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/filibsrc
    TMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/filibsrc-tmp
    STAMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/filibsrc-stamp
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/filibsrc
    BUILD_IN_SOURCE TRUE
    PATCH_COMMAND cp -r ${CMAKE_CURRENT_SOURCE_DIR}/filibsrc ${CMAKE_CURRENT_BINARY_DIR}
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/filibsrc/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/filibsrc-install
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/filibsrc-install/include)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/filibsrc-install/include DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/filibsrc-install/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/filibsrc-install/bin DESTINATION ${CMAKE_INSTALL_PREFIX})

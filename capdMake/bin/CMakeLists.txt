get_target_property(CAPD_CFLAGS ${PROJECT_NAME} COMPILE_OPTIONS)
get_target_property(CAPD_LIBS ${PROJECT_NAME} INTERFACE_LINK_LIBRARIES)
get_target_property(FILIB_CFLAGS filib COMPILE_OPTIONS)

string(REPLACE ";" " " CAPD_CFLAGS "${CAPD_CFLAGS}")
string(REPLACE ";" " " FILIB_CFLAGS "${FILIB_CFLAGS}")
string(JOIN " " CAPD_CFLAGS "-std=c++${CMAKE_CXX_STANDARD} ${CAPD_CFLAGS} ${FILIB_CFLAGS}")

string(REPLACE "filib" "-lfilib" CAPD_LIBS "${CAPD_LIBS}")
string(REPLACE ";" " " CAPD_LIBS "${CAPD_LIBS}")

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd-config.in
    ${CMAKE_BINARY_DIR}/bin/capd-config
    )

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd.pc.in
    ${CMAKE_BINARY_DIR}/bin/capd.pc
    )


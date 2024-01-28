get_target_property(CAPD_CFLAGS ${PROJECT_NAME} COMPILE_OPTIONS)
get_target_property(CAPD_LIBS ${PROJECT_NAME} INTERFACE_LINK_LIBRARIES)
get_target_property(CXXSTD ${PROJECT_NAME} CXX_STANDARD)

string(REPLACE ";" " " CAPD_CFLAGS "${CAPD_CFLAGS}")
string(JOIN " " CAPD_CFLAGS "-std=c++${CXXSTD}" "${CAPD_CFLAGS}")
string(REPLACE ";" " " CAPD_LIBS "${CAPD_LIBS}")
string(REPLACE "-lprim" "${CMAKE_BINARY_DIR}/capdExt/filibsrc/libprim/.libs/libprim.a" CAPD_LIBS "${CAPD_LIBS}")

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd-config.in
    ${CMAKE_BINARY_DIR}/bin/capd-config
    )

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd.pc.in
    ${CMAKE_BINARY_DIR}/bin/capd.pc
    )


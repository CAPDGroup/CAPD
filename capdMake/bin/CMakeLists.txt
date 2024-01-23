get_target_property(CAPD_CFLAGS ${PROJECT_NAME} COMPILE_OPTIONS)
get_target_property(CAPD_LIBS ${PROJECT_NAME} INTERFACE_LINK_LIBRARIES)

string(REPLACE ";" " " CAPD_CFLAGS "${CAPD_CFLAGS}")
string(REPLACE ";" " " CAPD_LIBS "${CAPD_LIBS}")

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd-config.in
    ${CMAKE_BINARY_DIR}/bin/capd-config
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ
    )

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd.pc.in
    ${CMAKE_BINARY_DIR}/bin/capd.pc
    FILE_PERMISSIONS OWNER_READ
    )


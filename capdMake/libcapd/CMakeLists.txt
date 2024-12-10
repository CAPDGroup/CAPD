###############################################################################################################################
# PROCESS CFLAGS
###############################################################################################################################

get_target_property(CAPD_CFLAGS ${PROJECT_NAME} COMPILE_OPTIONS)

string(REPLACE ";" " " CAPD_CFLAGS "${CAPD_CFLAGS}")
string(JOIN " " CAPD_CFLAGS "-std=c++${CMAKE_CXX_STANDARD} ${CAPD_CFLAGS}")

if(${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")
    
    get_target_property(FILIB_CFLAGS filib COMPILE_OPTIONS)
    string(REPLACE ";" " " FILIB_CFLAGS "${FILIB_CFLAGS}")
    string(JOIN " " CAPD_CFLAGS "${FILIB_CFLAGS}")

endif()

###############################################################################################################################
# PROCESS LIBS
###############################################################################################################################

get_target_property(CAPD_LIBS ${PROJECT_NAME} INTERFACE_LINK_LIBRARIES)

if(${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")

    string(REPLACE "filib" "-lfilib" CAPD_LIBS "${CAPD_LIBS}")

endif()

string(REPLACE ";" " " CAPD_LIBS "${CAPD_LIBS}")

###############################################################################################################################
# CONFIGURE capd-config 
# ( @todo fix duplication with another CMakeLists.txt )
###############################################################################################################################

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd-config.in
    ${CMAKE_CURRENT_BINARY_DIR}/capd-config
    @ONLY)

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/capd-config
    PERMISSIONS OWNER_EXECUTE OWNER_READ WORLD_READ WORLD_EXECUTE
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

###############################################################################################################################
# CONFIGURE capd.pc
###############################################################################################################################

set(PKG_CONFIG_PATH ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/capd.pc
    @ONLY)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/capd.pc
  PERMISSIONS OWNER_READ WORLD_READ
  DESTINATION ${PKG_CONFIG_PATH})

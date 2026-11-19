###############################################################################################################################
# Process cflags
###############################################################################################################################

get_target_property(CAPD_CFLAGS ${PROJECT_NAME} COMPILE_OPTIONS)

string(REPLACE ";" " " CAPD_CFLAGS "${CAPD_CFLAGS}")
string(JOIN " " CAPD_CFLAGS "-std=c++${CMAKE_CXX_STANDARD} ${CAPD_CFLAGS}")

if(${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")
    
    get_target_property(FILIB_CFLAGS filib COMPILE_OPTIONS)
    string(REPLACE ";" " " FILIB_CFLAGS "${FILIB_CFLAGS}")
    string(JOIN " " CAPD_CFLAGS "${CAPD_CFLAGS} ${FILIB_CFLAGS}")

    set(CAPD_EXT_INCLUDES "-I${PROJECT_SOURCE_DIR}/capdExt/filibsrc")

elseif(${CAPD_INTERVAL_TYPE} STREQUAL "TIGHT")

    # include directories and compiler flags of TIGHT and its dependencies (NFG, SIMDe, CORE-MATH)
    set(CAPD_EXT_INCLUDES_LIST "-I${CAPD_TIGHT_DIR}/include")
    foreach(_t nfg_iface coremath)
        get_target_property(_inc ${_t} INTERFACE_INCLUDE_DIRECTORIES)
        foreach(_d ${_inc})
            string(REGEX REPLACE "^\\$<BUILD_INTERFACE:(.*)>$" "\\1" _d "${_d}")
            list(APPEND CAPD_EXT_INCLUDES_LIST "-I${_d}")
        endforeach()
    endforeach()
    get_target_property(_nfgflags nfg_iface INTERFACE_COMPILE_OPTIONS)
    list(JOIN CAPD_EXT_INCLUDES_LIST " " CAPD_EXT_INCLUDES)
    list(JOIN _nfgflags " " _nfgflags)
    string(JOIN " " CAPD_CFLAGS "${CAPD_CFLAGS} ${_nfgflags}")

elseif(${CAPD_INTERVAL_TYPE} STREQUAL "CXSC")

    set(CAPD_EXT_INCLUDES_LIST
        "-I${PROJECT_SOURCE_DIR}/capdExt/cxsc-2-5-4/src"
        "-I${PROJECT_SOURCE_DIR}/capdExt/cxsc-2-5-4/src/fi_lib"
        "-I${PROJECT_SOURCE_DIR}/capdExt/cxsc-2-5-4/src/rts"
        "-I${PROJECT_SOURCE_DIR}/capdExt/cxsc-2-5-4/src/asm"
        "-I${CMAKE_BINARY_DIR}/capdExt/cxsc-2-5-4"
    )    

    string(REPLACE ";" " " CAPD_EXT_INCLUDES "${CAPD_EXT_INCLUDES_LIST}")

endif()

###############################################################################################################################
# Process libs
###############################################################################################################################

if(${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")

    get_target_property(CAPD_LIBS ${PROJECT_NAME} INTERFACE_LINK_LIBRARIES)
    string(REPLACE "filib" "-lfilib" CAPD_LIBS "${CAPD_LIBS}")
    string(REPLACE ";" " " CAPD_LIBS "${CAPD_LIBS}")
    set(CAPD_EXT_LIB_DIR "-L${CMAKE_BINARY_DIR}/capdExt/filibsrc")

elseif(${CAPD_INTERVAL_TYPE} STREQUAL "TIGHT")

    set(CAPD_LIBS "-ltight -lcoremath")
    set(CAPD_EXT_LIB_DIR "-L${CMAKE_BINARY_DIR}/tight -L${CMAKE_BINARY_DIR}/_deps/core-math-cpp-build")

elseif(${CAPD_INTERVAL_TYPE} STREQUAL "CXSC")

    get_target_property(CAPD_LIBS ${PROJECT_NAME} INTERFACE_LINK_LIBRARIES)
    string(REPLACE "cxsc" "-lcxsc" CAPD_LIBS "${CAPD_LIBS}")
    string(REPLACE ";" " " CAPD_LIBS "${CAPD_LIBS}")
    set(CAPD_EXT_LIB_DIR "-L${CMAKE_BINARY_DIR}/capdExt/cxsc-2-5-4")

endif()

###############################################################################################################################
# Configure (local)
###############################################################################################################################

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd-config.local.in
    ${CMAKE_BINARY_DIR}/bin/capd-config
    @ONLY)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd.pc.local.in
    ${CMAKE_BINARY_DIR}/bin/capd.pc
    @ONLY)

###############################################################################################################################
# Configure (installation)
###############################################################################################################################

set(PKG_CONFIG_PATH ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd-config.in
    ${CMAKE_CURRENT_BINARY_DIR}/capd-config
    @ONLY)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/capd.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/capd.pc
    @ONLY)

###############################################################################################################################
# Install
###############################################################################################################################

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/capd-config
    PERMISSIONS OWNER_EXECUTE OWNER_READ WORLD_READ WORLD_EXECUTE
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/capd.pc
  PERMISSIONS OWNER_READ WORLD_READ
  DESTINATION ${PKG_CONFIG_PATH})

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

add_subdirectory(include)

add_subdirectory(src/capd/alglib)

if (${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")

    add_subdirectory(filibsrc)

elseif(${CAPD_INTERVAL_TYPE} STREQUAL "CXSC")

    set(BUILD_SHARED OFF CACHE BOOL "Build the shared version of the library" FORCE)
    add_subdirectory(cxsc-2-5-4)

endif()



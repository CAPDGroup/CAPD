target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

add_subdirectory(include)

add_subdirectory(src/capd/alglib)

if(${CAPD_INTERVAL_TYPE} STREQUAL "FILIB")

    add_subdirectory(filibsrc)
    
endif()


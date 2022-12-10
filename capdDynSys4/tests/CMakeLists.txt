# -----------------------------------------------------------------------------------------------------------------------------
# Standalone test executables
# -----------------------------------------------------------------------------------------------------------------------------
add_subdirectory(nonAutTest)
add_subdirectory(BasicCnTaylorTest)
add_subdirectory(diffIncl)

add_subdirectory(BasicCurveMaskTest)

if (false)

    set (TEST_PROJECT_NAME capdDynSys4_Test)
    message("${TEST_PROJECT_NAME}")

    add_executable(${TEST_PROJECT_NAME})

    target_link_directories(${TEST_PROJECT_NAME} PRIVATE /media/avert/Pergamon/whiteeye/CAPD_Experimental/build)
    target_link_libraries(${TEST_PROJECT_NAME} capd_experimental)

    add_subdirectory(BasicCurveMaskTest)




    find_package (Boost REQUIRED COMPONENTS unit_test_framework)
    # target_link_directories(${TEST_PROJECT_NAME} PRIVATE ${Boost_LIBRARY_DIRS})
    target_link_libraries (${TEST_PROJECT_NAME} Boost::unit_test_framework)

endif()

build_test_executable_monofile(matrixMDTest)
build_test_executable_monofile(vectMDTest)
build_test_executable_monofile(vecttst)

build_boost_test_executable(vectalgTest vectalgTestMain.cpp StaticArrayInitTest.cpp)
build_boost_test_executable(C11vectalgTest C11vectalgTestMain.cpp C11vectalgInitTest.cpp )


# target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/vectalgTestMain.cpp)
# target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/vectMDTest.cpp)
# target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/C11vectalgInitTest.cpp)
# target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/vecttst.cpp)
# target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/matrixMDTest.cpp)
# target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StaticArrayInitTest.cpp)
# target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/C11vectalgTestMain.cpp)

function(build_test_executable TEST_PROJECT_NAME)

  add_executable(${TEST_PROJECT_NAME})
  target_compile_options(${TEST_PROJECT_NAME} PRIVATE -O0 -DFILIB_EXTENDED)

  target_link_libraries (${TEST_PROJECT_NAME} filib)
  target_sources(${TEST_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PROJECT_NAME}.cpp)

  add_test(NAME ${TEST_PROJECT_NAME} COMMAND ${TEST_PROJECT_NAME})

endfunction()

build_test_executable(check_0f) # native_directed
build_test_executable(check_0t) # native_switched
# build_test_executable(check_1f) # multiplicative
# build_test_executable(check_2f) # no_rounding
# build_test_executable(check_4f) # native_onesided_global
build_test_executable(check_5f) # pred_succ_rounding



##
# HALMD test tools
#
add_library(test_tools_lua
  lua.cpp
)

add_library(test_tools_ctest
  ctest.cpp
)

add_executable(test_tools_ctest_linking
  ctest_linking.cpp
)
target_link_libraries(test_tools_ctest_linking
  ${HALMD_TEST_LIBRARIES}
)

set(TEST_TOOLS_CTEST_LINKING_LOCATION_FILE
  ${CMAKE_CURRENT_BINARY_DIR}/test_tools_ctest_linking_location)
file(GENERATE OUTPUT ${TEST_TOOLS_CTEST_LINKING_LOCATION_FILE}
  CONTENT "$<TARGET_FILE:test_tools_ctest_linking>")

configure_file(
  ctest_linking.cmake.in
  test_tools_ctest_linking.cmake
  @ONLY
)
add_test(tools/ctest_linking
  test_tools_ctest_linking.cmake
)

add_executable(test_tools_init
  init.cpp
  init_duplicate.cpp
)
target_link_libraries(test_tools_init
  ${HALMD_TEST_LIBRARIES}
)
add_test(tools/init
  test_tools_init --log_level=test_suite
)

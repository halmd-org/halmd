add_library(test_tools_lua
  lua.cpp
)
target_link_libraries(test_tools_lua
  ${LUABIND_LIBRARY}
  ${LUA_LIBRARIES}
  dl
)

# The test_tools_ctest library contains a BOOST_GLOBAL_FIXTURE, which is not
# called explicitly by the executable it is linked against. To prevent the
# linker from optimizing away the fixture, we create a shared library.
add_library(test_tools_ctest SHARED
  ctest.cpp
)
target_link_libraries(test_tools_ctest
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
# The jfs2 filesystem on AIX has the nerving property that
# one cannot overwrite a shared library that has been used
# previously, even if the executable is not running anymore.
# We therefore delete the target output before linking.
get_target_property(
  TEST_TOOLS_CTEST_LOCATION test_tools_ctest LOCATION
)
add_custom_command(TARGET test_tools_ctest PRE_LINK
  COMMAND "${CMAKE_COMMAND}" -E remove -f "${TEST_TOOLS_CTEST_LOCATION}"
)

add_executable(test_tools_ctest_linking
  ctest_linking.cpp
)
target_link_libraries(test_tools_ctest_linking
  ${HALMD_TEST_LIBRARIES}
)
get_target_property(
  TEST_TOOLS_CTEST_LINKING_LOCATION test_tools_ctest_linking LOCATION
)
configure_file(
  ctest_linking.cmake.in
  test_tools_ctest_linking.cmake
  @ONLY
)
add_test(ctest_linking
  test_tools_ctest_linking.cmake
)

add_executable(test_init
  init.cpp
  init_duplicate.cpp
)
target_link_libraries(test_init
  ${HALMD_TEST_LIBRARIES}
)
add_test(init
  test_init --log_level=test_suite
)

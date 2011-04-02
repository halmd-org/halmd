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

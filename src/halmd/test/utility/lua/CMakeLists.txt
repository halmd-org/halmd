add_executable(test_lua_registry
  lua_registry.cpp
)
target_link_libraries(test_lua_registry
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(lua_registry
  test_lua_registry --log_level=test_suite
)

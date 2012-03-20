add_executable(test_unit_utility_lua_dependency_injection
  dependency_injection.cpp
)
target_link_libraries(test_unit_utility_lua_dependency_injection
  test_tools_lua
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/lua/dependency_injection
  test_unit_utility_lua_dependency_injection --log_level=test_suite
)

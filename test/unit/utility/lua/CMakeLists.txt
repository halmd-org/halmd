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
add_executable(test_unit_utility_lua_function
  function.cpp
)
target_link_libraries(test_unit_utility_lua_function
  test_tools_lua
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/lua/function
  test_unit_utility_lua_function --log_level=test_suite
)
add_executable(test_unit_utility_lua_map_converter
  map_converter.cpp
)
target_link_libraries(test_unit_utility_lua_map_converter
  halmd_utility_lua
  test_tools_lua
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/lua/map_converter
  test_unit_utility_lua_map_converter --log_level=test_suite
)

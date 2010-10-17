add_executable(test_dependency_injection
  dependency_injection.cpp
)
target_link_libraries(test_dependency_injection
  test_tools_lua
  ${HALMD_COMMON_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(dependency_injection
  test_dependency_injection --log_level=test_suite
)
add_executable(test_map_converter
  map_converter.cpp
)
target_link_libraries(test_map_converter
  test_tools_lua
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(map_converter
  test_map_converter --log_level=test_suite
)

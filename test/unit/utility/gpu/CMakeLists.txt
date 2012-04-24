add_executable(test_unit_utility_shared_memory
  shared_memory.cpp
)
target_link_libraries(test_unit_utility_shared_memory
  pthread
  rt
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(unit/utility/shared_memory
  test_unit_utility_shared_memory --log_level=message
)

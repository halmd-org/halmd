add_executable(test_typed_value
  typed_value.cpp
)
target_link_libraries(test_typed_value
  halmd_utility_program_options
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(typed_value
  test_typed_value --log_level=test_suite
)
add_executable(test_accumulating_value
  accumulating_value.cpp
)
target_link_libraries(test_accumulating_value
  halmd_utility_program_options
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(accumulating_value
  test_accumulating_value --log_level=test_suite
)

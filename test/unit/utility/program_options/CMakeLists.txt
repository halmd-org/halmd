add_executable(test_unit_utility_program_options_typed_value
  typed_value.cpp
)
target_link_libraries(test_unit_utility_program_options_typed_value
  halmd_utility_program_options
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/program_options/typed_value
  test_unit_utility_program_options_typed_value --log_level=test_suite
)
add_executable(test_unit_utility_program_options_accumulating_value
  accumulating_value.cpp
)
target_link_libraries(test_unit_utility_program_options_accumulating_value
  halmd_utility_program_options
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/program_options/accumulating_value
  test_unit_utility_program_options_accumulating_value --log_level=test_suite
)

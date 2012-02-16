add_executable(test_unit_numeric_fixed_vector
  fixed_vector.cpp
)
target_link_libraries(test_unit_numeric_fixed_vector
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/numeric/fixed_vector
  test_unit_numeric_fixed_vector --log_level=test_suite
)

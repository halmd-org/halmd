add_executable(test_unit_h5xx
  h5xx.cpp
)
target_link_libraries(test_unit_h5xx
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/h5xx/attribute
  test_unit_h5xx --run_test=h5xx_attribute --log_level=test_suite
)
add_test(unit/h5xx/dataset
  test_unit_h5xx --run_test=h5xx_dataset --log_level=test_suite
)
add_test(unit/h5xx/group
  test_unit_h5xx --run_test=h5xx_group --log_level=test_suite
)

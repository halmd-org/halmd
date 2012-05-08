add_executable(test_unit_algorithm_multi_range
  multi_range.cpp
)
target_link_libraries(test_unit_algorithm_multi_range
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/algorithm/multi_range
  test_unit_algorithm_multi_range --log_level=test_suite
)

if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif(HALMD_WITH_GPU)
add_subdirectory(host)

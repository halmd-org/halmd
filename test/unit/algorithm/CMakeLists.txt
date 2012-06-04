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

add_executable(test_unit_algorithm_permute
  permute.cpp
)
target_link_libraries(test_unit_algorithm_permute
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/algorithm/permute
  test_unit_algorithm_permute --log_level=test_suite
)

add_executable(test_unit_algorithm_radix_sort
  radix_sort.cpp
)
target_link_libraries(test_unit_algorithm_radix_sort
  ${HALMD_TEST_LIBRARIES}
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_algorithm_radix_sort
    halmd_algorithm_gpu
  )
endif()
add_test(unit/algorithm/radix_sort
  test_unit_algorithm_radix_sort --log_level=test_suite
)

if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif(HALMD_WITH_GPU)
add_subdirectory(host)

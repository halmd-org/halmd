add_executable(test_unit_algorithm_gpu_apply
  apply.cpp
  apply_kernel.cu
)
target_link_libraries(test_unit_algorithm_gpu_apply
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/algorithm/gpu/apply
  test_unit_algorithm_gpu_apply --log_level=test_suite
)

add_executable(test_unit_algorithm_gpu_reduce
  reduce.cpp
  reduce_kernel.cu
)
target_link_libraries(test_unit_algorithm_gpu_reduce
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/algorithm/gpu/reduce
  test_unit_algorithm_gpu_reduce --log_level=test_suite
)

add_executable(test_unit_algorithm_scan
  scan_kernel.cu
  scan.cpp
)
target_link_libraries(test_unit_algorithm_scan
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/algorithm/gpu/scan
  test_unit_algorithm_scan --log_level=test_suite
)

add_executable(test_unit_algorithm_gpu_copy_if
  copy_if_kernel.cu
  copy_if.cpp
)
target_link_libraries(test_unit_algorithm_gpu_copy_if
  ${HALMD_TEST_LIBRARIES}
  halmd_utility_gpu
)
add_test(unit/algorithm/gpu/copy_if
  test_unit_algorithm_gpu_copy_if --log_level=test_suite
)

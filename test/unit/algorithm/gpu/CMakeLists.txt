add_executable(test_unit_algorithm_gpu_apply
  apply.cpp
  apply_kernel.cu
)
target_link_libraries(test_unit_algorithm_gpu_apply
  ${HALMD_TEST_LIBRARIES}
  ${HALMD_COMMON_LIBRARIES}
  ${CUDA_LIBRARIES}
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
  ${Boost_SYSTEM_LIBRARY}
  ${CUDA_LIBRARIES}
)
add_test(unit/algorithm/gpu/scan
  test_unit_algorithm_scan --log_level=test_suite
)

add_executable(test_unit_algorithm_radix_sort
  radix_sort.cpp
)
target_link_libraries(test_unit_algorithm_radix_sort
  halmd_algorithm_gpu
  halmd_random_gpu
  ${HALMD_TEST_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${CUDA_LIBRARIES}
)
add_test(unit/algorithm/gpu/radix_sort
  test_unit_algorithm_radix_sort --log_level=test_suite
)

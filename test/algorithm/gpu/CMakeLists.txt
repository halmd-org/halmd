add_executable(test_gpu_reduce
  reduce.cpp
  reduce_kernel.cu
)
target_link_libraries(test_gpu_reduce
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  ${HALMD_COMMON_LIBRARIES}
  ${CUDA_LIBRARIES}
)
add_test(algorithm_gpu_reduce
  test_gpu_reduce --log_level=test_suite
)

add_executable(test_scan
  scan_kernel.cu
  scan.cpp
)
target_link_libraries(test_scan
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${CUDA_LIBRARIES}
)
add_test(algorithm_gpu_scan
  test_scan --log_level=test_suite
)

add_executable(test_radix_sort
  radix_sort.cpp
)
target_link_libraries(test_radix_sort
  halmd_algorithm_gpu
  halmd_random_gpu
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${CUDA_LIBRARIES}
)
add_test(algorithm_gpu_radix_sort
  test_radix_sort --log_level=test_suite
)

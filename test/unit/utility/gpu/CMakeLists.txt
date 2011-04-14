add_executable(test_unit_utility_gpu_variant
  variant.cpp
  variant_kernel.cu
)
target_link_libraries(test_unit_utility_gpu_variant
  ${HALMD_TEST_LIBRARIES}
  ${CUDA_LIBRARIES}
)
add_test(unit/utility/gpu/variant
  test_unit_utility_gpu_variant --log_level=test_suite
)


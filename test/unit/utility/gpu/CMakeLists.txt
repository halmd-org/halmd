# dsfloat_cuda_vector test
cuda_add_executable(test_unit_utility_gpu_dsfloat_cuda_vector
  dsfloat_cuda_vector.cpp
  dsfloat_cuda_vector_kernel.cu
)
target_link_libraries(test_unit_utility_gpu_dsfloat_cuda_vector
  ${HALMD_TEST_LIBRARIES}
)
halmd_add_gpu_test(unit/utility/gpu/dsfloat_cuda_vector
  test_unit_utility_gpu_dsfloat_cuda_vector --log_level=test_suite
)

# shared_memory test
add_executable(test_unit_utility_gpu_shared_memory
  shared_memory.cpp
)
target_link_libraries(test_unit_utility_gpu_shared_memory
  ${HALMD_TEST_LIBRARIES}
)
halmd_add_gpu_test(unit/utility/gpu/shared_memory
  test_unit_utility_gpu_shared_memory --log_level=test_suite
)

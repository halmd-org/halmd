add_executable(test_unit_numeric_fixed_vector
  fixed_vector.cpp
)
target_link_libraries(test_unit_numeric_fixed_vector
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/numeric/fixed_vector
  test_unit_numeric_fixed_vector --log_level=test_suite
)
if(HALMD_WITH_GPU)
  cuda_add_executable(test_unit_numeric_fixed_vector_cuda_vector_converter
    fixed_vector_cuda_vector_converter.cpp
    fixed_vector_cuda_vector_converter_kernel.cu
  )
  target_link_libraries(test_unit_numeric_fixed_vector_cuda_vector_converter
    ${HALMD_TEST_LIBRARIES}
    ${CUDA_LIBRARIES}
  )
  halmd_add_gpu_test(unit/numeric/fixed_vector/cuda_vector_converter
    test_unit_numeric_fixed_vector_cuda_vector_converter --log_level=test_suite
  )
endif(HALMD_WITH_GPU)

if(HALMD_WITH_GPU)
  # dsfun test
  add_executable(test_unit_numeric_mp_dsfun
    dsfun.cpp
    dsfun_kernel.cu
  )
  target_link_libraries(test_unit_numeric_mp_dsfun
    ${HALMD_TEST_LIBRARIES}
  )
  halmd_add_gpu_test(unit/numeric/mp/dsfun
    test_unit_numeric_mp_dsfun --log_level=test_suite
  )

  # dsfloat_cuda_vector test
  add_executable(test_unit_numeric_mp_dsfloat_cuda_vector
    dsfloat_cuda_vector.cpp
    dsfloat_cuda_vector_kernel.cu
  )
  target_link_libraries(test_unit_numeric_mp_dsfloat_cuda_vector
    ${HALMD_TEST_LIBRARIES}
  )
  halmd_add_gpu_test(unit/numeric/mp/dsfloat_cuda_vector
    test_unit_numeric_mp_dsfloat_cuda_vector --log_level=test_suite
  )
endif(HALMD_WITH_GPU)

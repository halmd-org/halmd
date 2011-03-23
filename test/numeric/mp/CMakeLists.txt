if(WITH_CUDA)
  add_executable(test_dsfun
    dsfun.cpp
    dsfun_kernel.cu
  )
  target_link_libraries(test_dsfun
    ${HALMD_TEST_LIBRARIES}
    ${CUDA_LIBRARIES}
  )
  add_test(numeric_mp_dsfun
    test_dsfun --log_level=test_suite
  )
endif(WITH_CUDA)

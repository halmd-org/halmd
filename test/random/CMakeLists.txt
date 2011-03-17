if(WITH_CUDA)
  add_subdirectory(gpu)

  add_executable(test_rand_dist rand_dist.cpp)
  target_link_libraries(test_rand_dist
    halmd_random_gpu
    ${HALMD_TEST_LIBRARIES}
    ${CUDA_LIBRARIES}
    ${GSL_LIBRARIES}
  )
  add_test(random_distributions
    test_rand_dist --log_level=test_suite
  )
endif(WITH_CUDA)

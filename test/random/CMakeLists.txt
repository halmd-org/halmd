add_executable(test_rand_dist
  rand_dist.cpp
)
if(WITH_CUDA)
  target_link_libraries(test_rand_dist
    halmd_random_gpu_random
  )
endif(WITH_CUDA)
target_link_libraries(test_rand_dist
  halmd_random_host_random
  halmd_io_logger
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(random_distributions
  test_rand_dist --log_level=test_suite
)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

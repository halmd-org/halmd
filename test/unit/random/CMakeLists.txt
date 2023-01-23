add_executable(test_unit_random_distributions
  distributions.cpp
)

if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_random_distributions
    halmd_random_gpu
  )
endif(HALMD_WITH_GPU)

target_link_libraries(test_unit_random_distributions
  halmd_random_host
  halmd_io
  ${HALMD_TEST_LIBRARIES}
)

add_test(unit/random/distributions/host/uniform
  test_unit_random_distributions --run_test=test_host_uniform\* --log_level=test_suite
)
add_test(unit/random/distributions/host/normal
  test_unit_random_distributions --run_test=test_host_normal\* --log_level=test_suite
)

if(HALMD_WITH_GPU)
  add_test(unit/random/distributions/gpu/uniform
    test_unit_random_distributions --run_test=test_gpu_uniform\* --log_level=test_suite
  )
  add_test(unit/random/distributions/gpu/normal
    test_unit_random_distributions --run_test=test_gpu_normal\* --log_level=test_suite
  )
endif(HALMD_WITH_GPU)

if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif(HALMD_WITH_GPU)

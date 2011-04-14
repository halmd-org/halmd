add_executable(test_unit_random_distributions
  distributions.cpp
)
if(WITH_CUDA)
  target_link_libraries(test_unit_random_distributions
    halmd_random_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_unit_random_distributions
  halmd_random_host
  halmd_io
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/random/distributions
  test_unit_random_distributions --log_level=test_suite
)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

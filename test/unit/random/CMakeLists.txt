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
add_test(unit/random/distributions
  test_unit_random_distributions --log_level=test_suite
)

if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif(HALMD_WITH_GPU)

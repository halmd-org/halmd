add_executable(test_position
  lattice.cpp
)
target_link_libraries(test_position
  halmd_mdsim_host
  halmd_mdsim_host_position
  halmd_observables_host
  halmd_observables_host_samples
  halmd_random_host
)
if(WITH_CUDA)
  target_link_libraries(test_position
    halmd_mdsim_gpu
    halmd_mdsim_gpu_position
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_position
  halmd_io
  halmd_mdsim
  halmd_random
  halmd_utility
  ${HALMD_COMMON_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(module_mdsim_position_lattice_host_2d
  test_position --run_test=lattice/host/2d --log_level=test_suite
)
add_test(module_mdsim_position_lattice_host_3d
  test_position --run_test=lattice/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(module_mdsim_position_lattice_gpu_2d
    test_position --run_test=lattice/gpu/2d --log_level=test_suite
  )
  add_test(module_mdsim_position_lattice_gpu_3d
    test_position --run_test=lattice/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)

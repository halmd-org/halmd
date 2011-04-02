add_executable(test_position_lattice
  lattice.cpp
)
target_link_libraries(test_position_lattice
  halmd_mdsim_host_particle
  halmd_mdsim_host_positions_lattice
  halmd_observables_host_phase_space
  halmd_observables_host_samples_phase_space
  halmd_random_host_random
)
if(WITH_CUDA)
  target_link_libraries(test_position_lattice
    halmd_mdsim_gpu_particle
    halmd_mdsim_gpu_positions_lattice
    halmd_observables_gpu_phase_space
    halmd_observables_gpu_samples_phase_space
    halmd_random_gpu_random
    halmd_utility_gpu_device
  )
endif(WITH_CUDA)
target_link_libraries(test_position_lattice
  halmd_io_logger
  halmd_mdsim_box
  halmd_mdsim_particle
  halmd_random_random
  halmd_utility_profiler
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(module_mdsim_position_lattice_host_2d
  test_position_lattice --run_test=lattice/host/2d --log_level=test_suite
)
add_test(module_mdsim_position_lattice_host_3d
  test_position_lattice --run_test=lattice/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(module_mdsim_position_lattice_gpu_2d
    test_position_lattice --run_test=lattice/gpu/2d --log_level=test_suite
  )
  add_test(module_mdsim_position_lattice_gpu_3d
    test_position_lattice --run_test=lattice/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)

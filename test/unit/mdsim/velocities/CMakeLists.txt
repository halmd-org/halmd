add_executable(test_unit_mdsim_velocities_boltzmann
  boltzmann.cpp
)
target_link_libraries(test_unit_mdsim_velocities_boltzmann
  halmd_mdsim_host_forces
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_observables_host_samples
  halmd_observables_host
  halmd_random_host
)
if(WITH_CUDA)
  target_link_libraries(test_unit_mdsim_velocities_boltzmann
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_unit_mdsim_velocities_boltzmann
  halmd_io
  halmd_mdsim
  halmd_observables
  halmd_random
  halmd_utility
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/velocities/boltzmann/host/2d
  test_unit_mdsim_velocities_boltzmann --run_test=boltzmann/host/2d --log_level=test_suite
)
add_test(unit/mdsim/velocities/boltzmann/host/3d
  test_unit_mdsim_velocities_boltzmann --run_test=boltzmann/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/velocities/boltzmann/gpu/2d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann/gpu/2d --log_level=test_suite
  )
  add_test(unit/mdsim/velocities/boltzmann/gpu/3d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)

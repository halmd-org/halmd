add_executable(test_velocity_boltzmann
  boltzmann.cpp
)
target_link_libraries(test_velocity_boltzmann
  halmd_mdsim_host_forces_zero
  halmd_mdsim_host_velocities_boltzmann
  halmd_mdsim_host_velocity
  halmd_mdsim_host_particle
  halmd_observables_host_thermodynamics
  halmd_random_host_random
)
if(WITH_CUDA)
  target_link_libraries(test_velocity_boltzmann
    halmd_mdsim_gpu_forces_zero
    halmd_mdsim_gpu_velocities_boltzmann
    halmd_mdsim_gpu_velocity
    halmd_mdsim_gpu_particle
    halmd_observables_gpu_thermodynamics
    halmd_random_gpu_random
    halmd_utility_gpu_device
  )
endif(WITH_CUDA)
target_link_libraries(test_velocity_boltzmann
  halmd_io_logger
  halmd_mdsim_box
  halmd_mdsim_particle
  halmd_observables_thermodynamics
  halmd_random_random
  halmd_utility_profiler
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(module_mdsim_velocity_boltzmann_host_2d
  test_velocity_boltzmann --run_test=boltzmann/host/2d --log_level=test_suite
)
add_test(module_mdsim_velocity_boltzmann_host_3d
  test_velocity_boltzmann --run_test=boltzmann/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(module_mdsim_velocity_boltzmann_gpu_2d
    test_velocity_boltzmann --run_test=boltzmann/gpu/2d --log_level=test_suite
  )
  add_test(module_mdsim_velocity_boltzmann_gpu_3d
    test_velocity_boltzmann --run_test=boltzmann/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)

add_executable(test_unit_mdsim_velocities_boltzmann
  boltzmann.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_velocities_boltzmann
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_velocities_boltzmann
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_observables
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/velocities/boltzmann/host/2d
  test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_host_2d --log_level=test_suite
)
add_test(unit/mdsim/velocities/boltzmann/host/3d
  test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  halmd_add_gpu_test(unit/mdsim/velocities/boltzmann/gpu/float/2d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_gpu_float_2d --log_level=test_suite
  )
  halmd_add_gpu_test(unit/mdsim/velocities/boltzmann/gpu/dsfloat/2d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_gpu_dsfloat_2d --log_level=test_suite
  )
  halmd_add_gpu_test(unit/mdsim/velocities/boltzmann/gpu/float/3d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_gpu_float_3d --log_level=test_suite
  )
  halmd_add_gpu_test(unit/mdsim/velocities/boltzmann/gpu/dsfloat/3d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_gpu_dsfloat_3d --log_level=test_suite
  )
endif()

# thermodynamics module
add_executable(test_unit_observables_thermodynamics
  thermodynamics.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_observables_thermodynamics
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_neighbours
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_potentials
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_observables_thermodynamics
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_neighbours
  halmd_mdsim_host_positions
  halmd_mdsim_host_potentials
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_observables
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/observables/thermodynamics/host/2d
  test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_host_2d --log_level=test_suite
)
add_test(unit/observables/thermodynamics/host/3d
  test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/observables/thermodynamics/gpu/2d
    test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_gpu_2d --log_level=test_suite
  )
  add_test(unit/observables/thermodynamics/gpu/3d
    test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_gpu_3d --log_level=test_suite
  )
endif()

# ssf module
add_executable(test_unit_observables_ssf
  ssf.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_observables_ssf
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_observables_ssf
  halmd_mdsim_host_positions
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_observables_utility
  halmd_observables
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/observables/ssf/host/2d
  test_unit_observables_ssf --run_test=ssf_host_2d --log_level=test_suite
)
add_test(unit/observables/ssf/host/3d
  test_unit_observables_ssf --run_test=ssf_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/observables/ssf/gpu/2d
    test_unit_observables_ssf --run_test=ssf_gpu_2d --log_level=test_suite
  )
  add_test(unit/observables/ssf/gpu/3d
    test_unit_observables_ssf --run_test=ssf_gpu_3d --log_level=test_suite
  )
endif()

# phase space sampler
add_executable(test_unit_observables_phase_space
  phase_space.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_observables_phase_space
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_observables_phase_space
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/observables/phase_space/host/2d
  test_unit_observables_phase_space --run_test=phase_space_host_2d --log_level=test_suite
)
add_test(unit/observables/phase_space/host/3d
  test_unit_observables_phase_space --run_test=phase_space_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/observables/phase_space/gpu/host/2d
    test_unit_observables_phase_space --run_test=phase_space_gpu_host_2d --log_level=test_suite
  )
  add_test(unit/observables/phase_space/gpu/host/3d
    test_unit_observables_phase_space --run_test=phase_space_gpu_host_3d --log_level=test_suite
  )
  add_test(unit/observables/phase_space/gpu/gpu/2d
    test_unit_observables_phase_space --run_test=phase_space_gpu_gpu_2d --log_level=test_suite
  )
  add_test(unit/observables/phase_space/gpu/gpu/3d
    test_unit_observables_phase_space --run_test=phase_space_gpu_gpu_3d --log_level=test_suite
  )
endif()

add_executable(test_thermodynamics
  thermodynamics.cpp
)
target_link_libraries(test_thermodynamics
  halmd_mdsim_host_forces_lennard_jones
  halmd_mdsim_host_integrators_verlet
  halmd_mdsim_host_integrators_verlet_nvt_andersen
  halmd_mdsim_host_neighbour
  halmd_mdsim_host_positions_lattice
  halmd_mdsim_host_velocities_boltzmann
  halmd_mdsim_host_velocity
  halmd_mdsim_host_particle
  halmd_observables_host_thermodynamics
  halmd_random_host_random
)
if(WITH_CUDA)
  target_link_libraries(test_thermodynamics
    halmd_mdsim_gpu_forces_lennard_jones
    halmd_mdsim_gpu_forces_pair_trunc_potentials
    halmd_mdsim_gpu_integrators_verlet
    halmd_mdsim_gpu_integrators_verlet_nvt_andersen
    halmd_mdsim_gpu_neighbour
    halmd_mdsim_gpu_positions_lattice
    halmd_mdsim_gpu_velocities_boltzmann
    halmd_mdsim_gpu_velocity
    halmd_mdsim_gpu_particle
    halmd_observables_gpu_thermodynamics
    halmd_random_gpu_random
    halmd_utility_gpu_device
  )
endif(WITH_CUDA)
target_link_libraries(test_thermodynamics
  halmd_io_logger
  halmd_mdsim_box
  halmd_mdsim_core
  halmd_mdsim_particle
  halmd_observables_thermodynamics
  halmd_utility_profiler
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(thermodynamics_host_2d
  test_thermodynamics --run_test=host/2d/* --log_level=test_suite
)
add_test(thermodynamics_host_3d
  test_thermodynamics --run_test=host/3d/* --log_level=test_suite
)
if(WITH_CUDA)
  add_test(thermodynamics_gpu_2d
    test_thermodynamics --run_test=gpu/2d/* --log_level=test_suite
  )
  add_test(thermodynamics_gpu_3d
    test_thermodynamics --run_test=gpu/3d/* --log_level=test_suite
  )
endif(WITH_CUDA)

# ssf module
add_executable(test_observables_ssf
  ssf.cpp
)
target_link_libraries(test_observables_ssf
  halmd_mdsim_host_positions_lattice
  halmd_mdsim_host_particle
  halmd_observables_host_density_mode
  halmd_observables_host_phase_space
  halmd_observables_host_samples_phase_space
  halmd_random_host_random
)
if(WITH_CUDA)
  target_link_libraries(test_observables_ssf
    halmd_mdsim_gpu_positions_lattice
    halmd_mdsim_gpu_particle
    halmd_observables_gpu_density_mode
    halmd_observables_gpu_phase_space
    halmd_observables_gpu_samples_phase_space
    halmd_random_gpu_random
    halmd_utility_gpu_device
  )
endif(WITH_CUDA)
target_link_libraries(test_observables_ssf
  halmd_io_logger
  halmd_mdsim_box
  halmd_mdsim_particle
  halmd_observables_density_mode
  halmd_observables_phase_space
  halmd_observables_ssf
  halmd_observables_samples_density_mode
  halmd_observables_utility_wavevector
  halmd_utility_profiler
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(module_observables_ssf_host_2d
  test_observables_ssf --run_test=ssf/host/2d --log_level=test_suite
)
add_test(module_observables_ssf_host_3d
  test_observables_ssf --run_test=ssf/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(module_observables_ssf_gpu_2d
    test_observables_ssf --run_test=ssf/gpu/2d --log_level=test_suite
  )
  add_test(module_observables_ssf_gpu_3d
    test_observables_ssf --run_test=ssf/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)

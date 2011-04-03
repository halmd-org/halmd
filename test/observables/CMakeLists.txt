# thermodynamics module
add_executable(test_thermodynamics
  thermodynamics.cpp
)
target_link_libraries(test_thermodynamics
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_observables_host
  halmd_random_host
)
if(WITH_CUDA)
  target_link_libraries(test_thermodynamics
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_thermodynamics
  halmd_io
  halmd_mdsim
  halmd_observables
  halmd_random
  halmd_utility
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
  halmd_mdsim_host_positions
  halmd_mdsim_host
  halmd_observables_host
  halmd_observables_host_samples
  halmd_random_host
)
if(WITH_CUDA)
  target_link_libraries(test_observables_ssf
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_observables_gpu_samples
    halmd_random_gpu
    halmd_utility_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_observables_ssf
  halmd_io
  halmd_mdsim
  halmd_observables
  halmd_observables_samples
  halmd_observables_utility
  halmd_random
  halmd_utility
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

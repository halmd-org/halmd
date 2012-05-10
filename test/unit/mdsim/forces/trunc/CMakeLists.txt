add_executable(test_unit_mdsim_forces_trunc_local_r4
  local_r4.cpp
)
if(WITH_CUDA)
  target_link_libraries(test_unit_mdsim_forces_trunc_local_r4
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_potentials
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_forces_trunc_local_r4
  halmd_mdsim_host_potentials
  halmd_mdsim_forces_trunc
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/forces/trunc/local_r4/host
  test_unit_mdsim_forces_trunc_local_r4 --run_test=host/local_r4 --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/forces/trunc/local_r4/gpu
    test_unit_mdsim_forces_trunc_local_r4 --run_test=gpu/local_r4 --log_level=test_suite
  )
endif()

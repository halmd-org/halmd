add_executable(test_unit_mdsim_positions_lattice
  lattice.cpp
)
target_link_libraries(test_unit_mdsim_positions_lattice
  halmd_mdsim_host
  halmd_mdsim_host_positions
  halmd_observables_host
  halmd_observables_host_samples
  halmd_random_host
)
if(WITH_CUDA)
  target_link_libraries(test_unit_mdsim_positions_lattice
    halmd_mdsim_gpu
    halmd_mdsim_gpu_positions
    halmd_observables_gpu
    halmd_observables_gpu_samples
    halmd_random_gpu
    halmd_utility_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_unit_mdsim_positions_lattice
  halmd_io
  halmd_mdsim
  halmd_random
  halmd_utility
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/positions/lattice/host/2d
  test_unit_mdsim_positions_lattice --run_test=lattice/host/2d --log_level=test_suite
)
add_test(unit/mdsim/positions/lattice/host/3d
  test_unit_mdsim_positions_lattice --run_test=lattice/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/positions/lattice/gpu/2d
    test_unit_mdsim_positions_lattice --run_test=lattice/gpu/2d --log_level=test_suite
  )
  add_test(unit/mdsim/positions/lattice/gpu/3d
    test_unit_mdsim_positions_lattice --run_test=lattice/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)

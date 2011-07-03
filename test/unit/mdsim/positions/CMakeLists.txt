add_executable(test_unit_mdsim_positions_lattice
  lattice.cpp
)
target_link_libraries(test_unit_mdsim_positions_lattice
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/positions/lattice/host/2d
  test_unit_mdsim_positions_lattice --run_test=lattice_host_2d --log_level=test_suite
)
add_test(unit/mdsim/positions/lattice/host/3d
  test_unit_mdsim_positions_lattice --run_test=lattice_host_3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/positions/lattice/gpu/2d
    test_unit_mdsim_positions_lattice --run_test=lattice_gpu_2d --log_level=test_suite
  )
  add_test(unit/mdsim/positions/lattice/gpu/3d
    test_unit_mdsim_positions_lattice --run_test=lattice_gpu_3d --log_level=test_suite
  )
endif(WITH_CUDA)

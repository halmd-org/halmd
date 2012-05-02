add_executable(test_unit_mdsim_positions_lattice
  lattice.cpp
)
if(WITH_CUDA)
  target_link_libraries(test_unit_mdsim_positions_lattice
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_positions_lattice
  halmd_mdsim_host_positions
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_random_host
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
endif()

add_executable(test_unit_mdsim_positions_lattice_primitive
  lattice_primitive.cpp
)
target_link_libraries(test_unit_mdsim_positions_lattice_primitive
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/positions/lattice_primitive
  test_unit_mdsim_positions_lattice_primitive --log_level=test_suite
)

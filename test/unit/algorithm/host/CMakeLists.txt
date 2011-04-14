add_executable(test_unit_algorithm_host_pick_lattice_points
  pick_lattice_points.cpp
)
target_link_libraries(test_unit_algorithm_host_pick_lattice_points
  halmd_io
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/algorithm/host/pick_lattice_points
  test_unit_algorithm_host_pick_lattice_points --log_level=test_suite
)
set_property(TEST unit/algorithm/host/pick_lattice_points
  PROPERTY TIMEOUT 60
)

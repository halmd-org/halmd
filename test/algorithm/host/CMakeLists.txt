add_executable(test_pick_lattice_points
  pick_lattice_points.cpp
)
target_link_libraries(test_pick_lattice_points
  halmd_io_logger
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(observables_utility_pick_lattice_points
  test_pick_lattice_points --log_level=test_suite
)

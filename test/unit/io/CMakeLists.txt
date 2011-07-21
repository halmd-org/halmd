add_executable(test_unit_io_logger
  logger.cpp
)
target_link_libraries(test_unit_io_logger
  halmd_io
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/logger
  test_unit_io_logger --log_level=test_suite
)

add_executable(test_unit_io_trajectory_h5md
  trajectory_h5md.cpp
)
target_link_libraries(test_unit_io_trajectory_h5md
  halmd_io_trajectory_readers
  halmd_io_trajectory_writers
  halmd_io_trajectory
  halmd_mdsim
  halmd_observables_host_samples
  halmd_utility
  halmd_io
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/trajectory_h5md/2d
  test_unit_io_trajectory_h5md --run_test=trajectory_h5md/2d --log_level=test_suite
)
add_test(unit/io/trajectory_h5md/3d
  test_unit_io_trajectory_h5md --run_test=trajectory_h5md/3d --log_level=test_suite
)

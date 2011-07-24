add_executable(test_unit_io_h5md_file
  file.cpp
)
target_link_libraries(test_unit_io_h5md_file
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/h5md/file
  test_unit_io_h5md_file --log_level=test_suite
)

add_executable(test_unit_io_h5md_trajectory
  trajectory.cpp
)
target_link_libraries(test_unit_io_h5md_trajectory
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/h5md/trajectory/2d
  test_unit_io_h5md_trajectory --run_test=2d --log_level=test_suite
)
add_test(unit/io/h5md/trajectory/3d
  test_unit_io_h5md_trajectory --run_test=3d --log_level=test_suite
)

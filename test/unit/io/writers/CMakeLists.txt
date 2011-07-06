add_executable(test_unit_io_writers_h5md
  h5md.cpp
)
target_link_libraries(test_unit_io_writers_h5md
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/writers/h5md
  test_unit_io_writers_h5md --log_level=test_suite
)


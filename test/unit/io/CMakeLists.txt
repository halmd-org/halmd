add_executable(test_unit_io_logger
  logger.cpp
)
target_link_libraries(test_unit_io_logger
  halmd_io
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/logger
  test_unit_io_logger --log_level=test_suite
)

add_subdirectory(h5md)

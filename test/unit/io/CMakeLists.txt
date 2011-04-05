add_executable(test_unit_io_logger
  logger.cpp
)
target_link_libraries(test_unit_io_logger
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/logger
  test_unit_io_logger --log_level=test_suite
)

add_executable(test_unit_io_profiling_writers
  profiling_writers.cpp
)
target_link_libraries(test_unit_io_profiling_writers
  halmd_io_profiling_writers
  halmd_io_profiling
  halmd_utility
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/io/profiling_writers
  test_unit_io_profiling_writers --log_level=test_suite
)

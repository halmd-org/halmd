add_executable(test_logger
  logger.cpp
)
target_link_libraries(test_logger
  halmd_io_logger
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(io_logger
  test_logger --log_level=test_suite
)

add_executable(test_profiling_writers
  profiling_writers.cpp
)
target_link_libraries(test_profiling_writers
  halmd_io_logger
  halmd_io_profiling_writers_h5md
  halmd_io_profiling_writers_log
  halmd_utility_profiler
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(io_profiling_writers
  test_profiling_writers --log_level=test_suite
)

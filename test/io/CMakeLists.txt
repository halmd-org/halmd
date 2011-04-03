add_executable(test_logger
  logger.cpp
)
target_link_libraries(test_logger
  halmd_io
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
  halmd_io_profiling_writers
  halmd_io_profiling
  halmd_utility
  halmd_io
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(io_profiling_writers
  test_profiling_writers --log_level=test_suite
)

add_executable(test_unit_cuda_wrapper_vector
  vector.cpp
)
target_link_libraries(test_unit_cuda_wrapper_vector
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/cuda_wrapper/vector
  test_unit_cuda_wrapper_vector --log_level=test_suite
)

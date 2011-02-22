add_executable(test_cuda_vector
  vector.cpp
)
target_link_libraries(test_cuda_vector
  ${HALMD_COMMON_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(cuda_vector
  test_cuda_vector --log_level=test_suite
)

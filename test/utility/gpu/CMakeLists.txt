add_executable(test_variant
  variant.cpp
  variant_kernel.cu
)
target_link_libraries(test_variant
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  ${CUDA_LIBRARIES}
)
add_test(gpu_variant
  test_variant --log_level=test_suite
)


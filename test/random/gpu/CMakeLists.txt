add_executable(test_rand48
  rand48.cpp
  )
target_link_libraries(test_rand48
  halmd_random_gpu
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${CUDA_LIBRARIES}
)
add_test(random_gpu_rand48
  test_rand48 --log_level=test_suite
)

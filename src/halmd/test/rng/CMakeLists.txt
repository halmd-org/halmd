if(WITH_CUDA)
  add_executable(rand48_test rand48.cpp)
  target_link_libraries(rand48_test
    gpu_rand48
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )

  add_executable(test_rand_dist rand_dist.cpp)
  target_link_libraries(test_rand_dist
    gpu_rand48
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${CUDA_LIBRARIES}
    ${GSL_LIBRARIES}
    )
  add_test(random_distributions test_rand_dist)
endif(WITH_CUDA)

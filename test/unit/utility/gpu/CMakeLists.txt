# shared_memory test
add_executable(test_unit_utility_gpu_shared_memory
  shared_memory.cpp
)
target_link_libraries(test_unit_utility_gpu_shared_memory
  ${HALMD_TEST_LIBRARIES}
)
halmd_add_gpu_test(unit/utility/gpu/shared_memory
  test_unit_utility_gpu_shared_memory --log_level=test_suite
)

##
# HALMD performance tests
#
if(HALMD_WITH_GPU)
  add_executable(test_performance_function_calls
    function_calls.cpp
    function_calls_extern.cpp
    function_calls_extern_kernel.cu
  )
else(HALMD_WITH_GPU)
  add_executable(test_performance_function_calls
    function_calls.cpp
    function_calls_extern.cpp
  )
endif(HALMD_WITH_GPU)
target_link_libraries(test_performance_function_calls
  test_tools_lua
  ${HALMD_TEST_LIBRARIES}
)
add_test(performance/function_calls
  test_performance_function_calls --log_level=message
)

if(HALMD_WITH_GPU)
  add_executable(test_performance_memory_allocation
    memory_allocation.cpp
  )
  target_link_libraries(test_performance_memory_allocation
    ${HALMD_TEST_LIBRARIES}
  )
  halmd_add_gpu_test(performance/memory_allocation
    test_performance_memory_allocation --log_level=message
  )

  add_executable(test_performance_reduction
    reduction.cpp
    reduction_kernel.cu
  )
  target_link_libraries(test_performance_reduction
    ${HALMD_TEST_LIBRARIES}
  )
  halmd_add_gpu_test(performance/reduction
    test_performance_reduction --log_level=message
  )
endif(HALMD_WITH_GPU)

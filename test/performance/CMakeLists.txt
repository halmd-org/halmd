##
# HALMD performance tests
#
if(WITH_CUDA)
  add_executable(test_performance_function_calls
    function_calls.cpp
    function_calls_extern.cpp
    function_calls_extern_kernel.cu
  )
else(WITH_CUDA)
  add_executable(test_performance_function_calls
    function_calls.cpp
    function_calls_extern.cpp
    function_calls_extern_kernel.cu
  )
endif(WITH_CUDA)
target_link_libraries(test_performance_function_calls
  test_tools_lua
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(performance/function_calls
  test_performance_function_calls --log_level=message
)

if(WITH_CUDA)
  add_executable(test_performance_fill_vs_memset
    fill_vs_memset.cpp
    fill_kernel.cu
  )
else(WITH_CUDA)
  add_executable(test_performance_fill_vs_memset
    fill_vs_memset.cpp
  )
endif(WITH_CUDA)
target_link_libraries(test_performance_fill_vs_memset
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(performance/fill_vs_memset/host
  test_performance_fill_vs_memset --run_test=host --log_level=message
)
if(WITH_CUDA)
  add_test(performance/fill_vs_memset/gpu
    test_performance_fill_vs_memset --run_test=gpu --log_level=message
  )
endif(WITH_CUDA)

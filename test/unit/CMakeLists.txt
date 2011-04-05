##
# HALMD unit tests
#

# Link unit tests against the test_tools_logger library,
# which enables logging to console for HALMD modules.
set(HALMD_TEST_LIBRARIES
  ${HALMD_TEST_LIBRARIES}
  test_tools_logger
)

add_subdirectory(algorithm)
add_subdirectory(h5xx)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(numeric)
add_subdirectory(observables)
add_subdirectory(random)
add_subdirectory(utility)
if(WITH_CUDA)
  add_subdirectory(cuda_wrapper)
endif(WITH_CUDA)

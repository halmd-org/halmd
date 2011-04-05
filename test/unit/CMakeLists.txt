##
# HALMD unit tests
#
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

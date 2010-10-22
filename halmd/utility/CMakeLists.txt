add_library(halmd_utility
  hostname.cpp
  profiler.cpp
)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

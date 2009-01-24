add_library(sample
  H5param.cpp
  correlation.cpp
  energy.cpp
  perf.cpp
  trajectory.cpp
  )

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

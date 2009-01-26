add_library(sample
  H5param.cpp
  correlation.cpp
  energy.cpp
  perf.cpp
  trajectory.cpp
  )
target_link_libraries(sample
  ${CUDA_LIBRARIES}
  gpu_tcf
  )

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

add_library(sample_common
  H5param.cpp
  energy.cpp
  perf.cpp
  trajectory.cpp
  )
add_library(sample_host
  correlation.cpp
  )
set_target_properties(sample_host PROPERTIES
  COMPILE_FLAGS "-UWITH_CUDA"
  )
if(WITH_CUDA)
  add_library(sample_gpu
    correlation.cpp
    )
  target_link_libraries(sample_gpu
    ${CUDA_LIBRARIES}
    gpu_tcf
    )
  add_subdirectory(gpu)
endif(WITH_CUDA)

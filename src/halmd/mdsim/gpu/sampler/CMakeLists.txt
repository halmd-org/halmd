add_library(halmd_mdsim_gpu_sample
  trajectory.cpp
  trajectory_wrapper.cu
  )
target_link_libraries(halmd_mdsim_gpu_sample
  ${CUDA_LIBRARIES}
  )

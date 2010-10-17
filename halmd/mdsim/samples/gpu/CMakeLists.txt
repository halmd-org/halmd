add_library(halmd_mdsim_samples_gpu
  trajectory.cpp
  )
target_link_libraries(halmd_mdsim_samples_gpu
  ${CUDA_LIBRARIES}
  )

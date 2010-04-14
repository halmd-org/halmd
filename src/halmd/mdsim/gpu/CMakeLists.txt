add_library(halmd_mdsim_gpu
  particle.cpp
  )
target_link_libraries(halmd_mdsim_gpu
  ${CUDA_LIBRARIES}
  )

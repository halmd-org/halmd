add_library(halmd_mdsim_gpu_forces
  lj_wrapper.cu
#lj.cpp
  )
target_link_libraries(halmd_mdsim_gpu_forces
  ${CUDA_LIBRARIES}
  )

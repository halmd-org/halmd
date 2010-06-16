add_library(halmd_mdsim_gpu_integrators
  verlet_kernel.cu
  verlet.cpp
  )
target_link_libraries(halmd_mdsim_gpu_integrators
  ${CUDA_LIBRARIES}
  )

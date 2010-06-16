add_library(halmd_mdsim_gpu_integrator
  verlet_kernel.cu
  verlet.cpp
  )
target_link_libraries(halmd_mdsim_gpu_integrator
  ${CUDA_LIBRARIES}
  )

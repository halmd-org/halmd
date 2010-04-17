add_library(halmd_mdsim_gpu
  particle.cpp
  )
target_link_libraries(halmd_mdsim_gpu
  ${CUDA_LIBRARIES}
  halmd_mdsim_gpu_forces
  halmd_mdsim_gpu_integrator
  )

add_subdirectory(forces)
add_subdirectory(integrator)

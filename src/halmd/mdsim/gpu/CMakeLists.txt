add_library(halmd_mdsim_gpu
  neighbour_wrapper.cu
  particle.cpp
  )
target_link_libraries(halmd_mdsim_gpu
  ${CUDA_LIBRARIES}
  )

add_subdirectory(forces)
add_subdirectory(integrator)
add_subdirectory(sample)
add_subdirectory(sort)

add_library(halmd_utility_gpu
  device.cpp
  )
target_link_libraries(halmd_utility_gpu
  ${CUDA_LIBRARIES}
  )

add_library(halmd_random_gpu
  rand48_kernel.cu
  random.cpp
)
target_link_libraries(halmd_random_gpu
  halmd_algorithm_gpu
)

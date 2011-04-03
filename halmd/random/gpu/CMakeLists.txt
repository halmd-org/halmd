halmd_add_module(halmd_random_gpu_random
  rand48_kernel.cu
  random_kernel.cu
  random.cpp
)
target_link_libraries(halmd_random_gpu_random
  halmd_algorithm_gpu
)

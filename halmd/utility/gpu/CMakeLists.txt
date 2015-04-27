halmd_add_library(halmd_utility_gpu
  device.cpp
  caching_array.cu
)
halmd_add_modules(
  libhalmd_utility_gpu_device
)

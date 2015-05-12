halmd_add_library(halmd_utility_gpu
  device.cpp
  device.cu
)
halmd_add_modules(
  libhalmd_utility_gpu_device
)

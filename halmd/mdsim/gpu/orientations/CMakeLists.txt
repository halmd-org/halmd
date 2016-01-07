halmd_add_library(halmd_mdsim_gpu_orientations
  uniform.cpp
  uniform_kernel.cu
)
halmd_add_modules(
  libhalmd_mdsim_gpu_orientations_uniform
)

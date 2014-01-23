halmd_add_library(halmd_mdsim_gpu_mobilities
  oseen.cpp
  oseen_kernel.cu
)
halmd_add_modules(
  libhalmd_mdsim_gpu_mobilities_oseen
)

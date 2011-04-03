halmd_add_library(halmd_mdsim_gpu_sorts
  hilbert.cpp
  hilbert_kernel.cu
)
halmd_add_modules(
  libhalmd_mdsim_gpu_sorts_hilbert
)

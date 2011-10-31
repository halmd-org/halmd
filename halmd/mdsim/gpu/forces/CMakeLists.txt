halmd_add_library(halmd_mdsim_gpu_forces
  pair_trunc.cpp
  zero.cpp
)
halmd_add_modules(
  libhalmd_mdsim_gpu_forces_pair_trunc
  libhalmd_mdsim_gpu_forces_zero
)

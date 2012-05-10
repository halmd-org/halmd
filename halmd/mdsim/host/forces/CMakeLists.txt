halmd_add_library(halmd_mdsim_host_forces
  pair_full.cpp
  pair_trunc.cpp
)
halmd_add_modules(
  libhalmd_mdsim_host_forces_pair_full
  libhalmd_mdsim_host_forces_pair_trunc
)

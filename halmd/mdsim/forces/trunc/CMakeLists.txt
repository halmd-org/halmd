halmd_add_library(halmd_mdsim_forces_trunc
  discontinuous.cpp
  local_r4.cpp
)
halmd_add_modules(
  libhalmd_mdsim_forces_trunc_discontinuous
  libhalmd_mdsim_forces_trunc_local_r4
)

halmd_add_library(halmd_mdsim_host_positions
  excluded_volume.cpp
  lattice.cpp
)
halmd_add_modules(
  libhalmd_mdsim_host_positions_excluded_volume
  libhalmd_mdsim_host_positions_lattice
)

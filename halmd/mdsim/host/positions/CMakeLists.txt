halmd_add_library(halmd_mdsim_host_positions
  lattice.cpp
  phase_space.cpp
)
halmd_add_modules(
  libhalmd_mdsim_host_positions_lattice
  libhalmd_mdsim_host_positions_phase_space
)

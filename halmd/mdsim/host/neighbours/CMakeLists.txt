halmd_add_library(halmd_mdsim_host_neighbours
  from_binning.cpp
  from_particle.cpp
)
halmd_add_modules(
  libhalmd_mdsim_host_neighbours_from_binning
  libhalmd_mdsim_host_neighbours_from_particle
)

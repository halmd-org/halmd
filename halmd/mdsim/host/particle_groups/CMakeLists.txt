halmd_add_modules(
  libhalmd_mdsim_host_particle_groups_from_range
  libhalmd_mdsim_host_particle_groups_all
)
halmd_add_library(halmd_mdsim_host_particle_groups
  from_range.cpp
  all.cpp
)

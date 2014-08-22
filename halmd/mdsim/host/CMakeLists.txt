halmd_add_modules(
  libhalmd_mdsim_host_binning
  libhalmd_mdsim_host_max_displacement
  libhalmd_mdsim_host_neighbour
  libhalmd_mdsim_host_particle
  libhalmd_mdsim_host_particle_group
  libhalmd_mdsim_host_region
)

add_subdirectory(integrators)
add_subdirectory(neighbours)
add_subdirectory(particle_groups)
add_subdirectory(positions)
add_subdirectory(potentials)
add_subdirectory(sorts)
add_subdirectory(velocities)

halmd_add_library(halmd_mdsim_host
  binning.cpp
  max_displacement.cpp
  neighbour.cpp
  particle.cpp
  particle_group.cpp
  region.cpp
)

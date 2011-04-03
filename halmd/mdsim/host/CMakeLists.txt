halmd_add_library(halmd_mdsim_host
  force.cpp
  neighbour.cpp
  particle.cpp
  velocity.cpp
)
halmd_add_modules(
  libhalmd_mdsim_host_force
  libhalmd_mdsim_host_neighbour
  libhalmd_mdsim_host_particle
  libhalmd_mdsim_host_velocity
)

add_subdirectory(forces)
add_subdirectory(integrators)
add_subdirectory(positions)
add_subdirectory(sorts)
add_subdirectory(velocities)

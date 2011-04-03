halmd_add_module(halmd_mdsim_host_force
  force.cpp
)
halmd_add_module(halmd_mdsim_host_neighbour
  neighbour.cpp
)
halmd_add_module(halmd_mdsim_host_particle
  particle.cpp
)
halmd_add_module(halmd_mdsim_host_velocity
  velocity.cpp
)

add_subdirectory(forces)
add_subdirectory(integrators)
add_subdirectory(positions)
add_subdirectory(sorts)
add_subdirectory(velocities)

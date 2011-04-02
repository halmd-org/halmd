add_library(halmd_mdsim_host_force
  force.cpp
)
add_library(halmd_mdsim_host_neighbour
  neighbour.cpp
)
add_library(halmd_mdsim_host_particle
  particle.cpp
)
add_library(halmd_mdsim_host_velocity
  velocity.cpp
)

add_subdirectory(forces)
add_subdirectory(integrators)
add_subdirectory(positions)
add_subdirectory(sorts)
add_subdirectory(velocities)

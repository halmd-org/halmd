add_library(halmd_mdsim_host
  force.cpp
  neighbour.cpp
  particle.cpp
  velocity.cpp
  )

add_subdirectory(forces)
add_subdirectory(integrator)
add_subdirectory(position)
add_subdirectory(sampler)
add_subdirectory(sort)
add_subdirectory(velocities)

add_library(halmd_mdsim_host
  box.cpp
  force.cpp
  neighbor.cpp
  particle.cpp
  random.cpp
  thermodynamics.cpp
  )

add_subdirectory(forces)
add_subdirectory(integrator)
add_subdirectory(position)
add_subdirectory(sample)
add_subdirectory(sort)
add_subdirectory(velocity)

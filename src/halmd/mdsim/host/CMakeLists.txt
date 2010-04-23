add_library(halmd_mdsim_host
  box.cpp
  neighbor.cpp
  particle.cpp
  random.cpp
  )

add_subdirectory(forces)
add_subdirectory(integrator)
add_subdirectory(position)
add_subdirectory(sample)
add_subdirectory(sort)
add_subdirectory(velocity)

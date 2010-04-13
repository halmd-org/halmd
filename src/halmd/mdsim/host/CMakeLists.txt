add_library(halmd_mdsim_host
  neighbor.cpp
  particle.cpp
  random.cpp
  )
target_link_libraries(halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrator
  halmd_mdsim_host_position
  halmd_mdsim_host_sort
  halmd_mdsim_host_velocity
  )

add_subdirectory(forces)
add_subdirectory(integrator)
add_subdirectory(position)
add_subdirectory(sort)
add_subdirectory(velocity)

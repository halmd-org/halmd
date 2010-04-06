add_library(halmd_mdsim_host
  neighbor.cpp
  particle.cpp
  random.cpp
  )
target_link_libraries(halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_position
  )

add_subdirectory(forces)
add_subdirectory(position)

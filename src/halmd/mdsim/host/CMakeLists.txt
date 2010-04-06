add_library(halmd_mdsim_host
  neighbor.cpp
  particle.cpp
  random.cpp
  )
target_link_libraries(halmd_mdsim_host
  halmd_mdsim_host_forces
  )

add_subdirectory(forces)

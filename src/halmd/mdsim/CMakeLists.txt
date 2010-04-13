add_library(halmd_mdsim
  box.cpp
  core.cpp
  force.cpp
  integrator.cpp
  neighbor.cpp
  particle.cpp
  position.cpp
  random.cpp
  sort.cpp
  velocity.cpp
  )
target_link_libraries(halmd_mdsim
  halmd_mdsim_host
  )

add_subdirectory(backend)
add_subdirectory(host)

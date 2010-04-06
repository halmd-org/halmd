add_library(halmd_mdsim
  box.cpp
  particle.cpp
  random.cpp
  )
target_link_libraries(halmd_mdsim
  halmd_mdsim_host
  )

add_subdirectory(backend)
add_subdirectory(host)

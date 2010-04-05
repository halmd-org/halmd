add_library(halmd_mdsim
  box.cpp
  )
target_link_libraries(halmd_mdsim
  halmd_mdsim_host_forces
  )

add_subdirectory(backend)
add_subdirectory(host)

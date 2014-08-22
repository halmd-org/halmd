halmd_add_library(halmd_mdsim
  box.cpp
  clock.cpp
  core.cpp
)
halmd_add_modules(
  libhalmd_mdsim_box
  libhalmd_mdsim_core
  libhalmd_mdsim_clock
)

add_subdirectory(forces)
add_subdirectory(geometries)

if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif(HALMD_WITH_GPU)
add_subdirectory(host)

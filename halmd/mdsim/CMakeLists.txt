halmd_add_library(halmd_mdsim
  box.cpp
  clock.cpp
  core.cpp
  neighbour.cpp
  sort.cpp
)
halmd_add_modules(
  libhalmd_mdsim_box
  libhalmd_mdsim_core
  libhalmd_mdsim_clock
  libhalmd_mdsim_neighbour
  libhalmd_mdsim_sort
)

add_subdirectory(forces)

if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif(HALMD_WITH_GPU)
add_subdirectory(host)

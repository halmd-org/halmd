halmd_add_library(halmd_io_trajectory
  reader.cpp
)
halmd_add_modules(
  libhalmd_io_trajectory_reader
)

add_subdirectory(readers)

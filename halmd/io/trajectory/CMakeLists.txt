halmd_add_library(halmd_io_trajectory
  reader.cpp
  writer.cpp
)
halmd_add_modules(
  libhalmd_io_trajectory_reader
  libhalmd_io_trajectory_writer
)

add_subdirectory(readers)
add_subdirectory(writers)

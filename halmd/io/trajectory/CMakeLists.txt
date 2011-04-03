halmd_add_module(halmd_io_trajectory_reader
  reader.cpp
)
halmd_add_module(halmd_io_trajectory_writer
  writer.cpp
)

add_subdirectory(readers)
add_subdirectory(writers)

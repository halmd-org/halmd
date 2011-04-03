add_library(halmd_io_trajectory_reader
  reader.cpp
)
add_library(halmd_io_trajectory_writer
  writer.cpp
)

add_subdirectory(readers)
add_subdirectory(writers)

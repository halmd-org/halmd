add_library(halmd_io_trajectory
  reader.cpp
  writer.cpp
  )

add_subdirectory(readers)
add_subdirectory(writers)

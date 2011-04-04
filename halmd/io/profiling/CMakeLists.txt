halmd_add_library(halmd_io_profiling
  writer.cpp
)
halmd_add_modules(
  libhalmd_io_profiling_writer
)

add_subdirectory(writers)

halmd_add_library(halmd_io_statevars
  writer.cpp
)
halmd_add_modules(
  libhalmd_io_statevars_writer
)

add_subdirectory(writers)

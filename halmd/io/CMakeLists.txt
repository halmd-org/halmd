halmd_add_library(halmd_io
  logger.cpp
)
halmd_add_modules(
  libhalmd_io_logger
)

add_subdirectory(profiling)
add_subdirectory(trajectory)
add_subdirectory(writers)

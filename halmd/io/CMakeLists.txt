halmd_add_library(halmd_io
  logger.cpp
)
halmd_add_modules(
  libhalmd_io_logger
)

add_subdirectory(readers)
add_subdirectory(writers)

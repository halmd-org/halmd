halmd_add_library(halmd_io_readers_h5md
  file.cpp
  truncate.cpp
)
halmd_add_modules(
  libhalmd_io_readers_h5md_file
  libhalmd_io_readers_h5md_truncate
)

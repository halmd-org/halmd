halmd_add_library(halmd_utility
  hostname.cpp
  options_parser.cpp
  posix_signal.cpp
  profiler.cpp
)
halmd_add_modules(
  libhalmd_utility_options_parser
  libhalmd_utility_profiler
)

add_subdirectory(lua)
add_subdirectory(program_options)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

add_library(halmd_utility
  hostname.cpp
  options_parser.cpp
  posix_signal.cpp
  profiler.cpp
)

add_subdirectory(lua_wrapper)
add_subdirectory(program_options)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

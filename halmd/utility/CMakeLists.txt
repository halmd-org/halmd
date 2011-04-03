add_library(halmd_utility_hostname
  hostname.cpp
)
add_library(halmd_utility_options_parser
  options_parser.cpp
)
add_library(halmd_utility_posix_signal
  posix_signal.cpp
)
add_library(halmd_utility_profiler
  profiler.cpp
)

add_subdirectory(lua)
add_subdirectory(program_options)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

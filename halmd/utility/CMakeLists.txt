halmd_add_library(halmd_utility
  hostname.cpp
  posix_signal.cpp
  profiler.cpp
  timer_service.cpp
  version.cpp
)
halmd_add_modules(
  libhalmd_utility_posix_signal
  libhalmd_utility_profiler
  libhalmd_utility_timer_service
  libhalmd_utility_version
)

add_subdirectory(lua)
if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif(HALMD_WITH_GPU)

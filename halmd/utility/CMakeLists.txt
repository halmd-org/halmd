halmd_add_library(halmd_utility
  hostname.cpp
  posix_signal.cpp
  profiler.cpp
  timer_service.cpp
)
halmd_add_modules(
  libhalmd_utility_posix_signal
  libhalmd_utility_profiler
  libhalmd_utility_timer_service
)

add_subdirectory(lua)
add_subdirectory(predicates)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

add_library(halmd_observables
  density_modes.cpp
  observable.cpp
  ssf.cpp
  thermodynamics.cpp
  trajectory.cpp
)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
add_subdirectory(host)
add_subdirectory(samples)
add_subdirectory(utility)

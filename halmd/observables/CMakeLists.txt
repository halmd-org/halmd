add_library(halmd_observables
  density_mode.cpp
  observable.cpp
  ssf.cpp
  thermodynamics.cpp
  phase_space.cpp
)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
add_subdirectory(host)
add_subdirectory(samples)
add_subdirectory(utility)

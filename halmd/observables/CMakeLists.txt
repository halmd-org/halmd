add_library(halmd_observables
  density_mode.cpp
  observable.cpp
  phase_space.cpp
  ssf.cpp
  thermodynamics.cpp
)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
add_subdirectory(host)
add_subdirectory(samples)
add_subdirectory(utility)

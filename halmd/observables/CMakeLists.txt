add_library(halmd_observables
  observable.cpp
  thermodynamics.cpp
  trajectory.cpp
)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
add_subdirectory(host)

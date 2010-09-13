add_library(halmd_observables
  thermodynamics.cpp
  )

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
add_subdirectory(host)

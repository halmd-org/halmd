add_library(halmd_observables_host
  mean_square_displacement.cpp
  thermodynamics.cpp
  trajectory.cpp
)
add_subdirectory(samples)

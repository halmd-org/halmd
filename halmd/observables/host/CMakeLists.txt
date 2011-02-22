add_library(halmd_observables_host
  mean_square_displacement.cpp
  thermodynamics.cpp
  trajectory.cpp
  velocity_autocorrelation.cpp
)
add_subdirectory(samples)

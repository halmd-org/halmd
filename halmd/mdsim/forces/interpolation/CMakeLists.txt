halmd_add_library(halmd_mdsim_forces_interpolation
  cubic_hermite.cpp
  linear.cpp
)
halmd_add_modules(
  libhalmd_mdsim_forces_interpolation_cubic_hermite
  libhalmd_mdsim_forces_interpolation_linear
)

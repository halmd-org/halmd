halmd_add_library(halmd_observables_gpu_fields
  density.cpp
  density_kernel.cu
)
halmd_add_modules(
  libhalmd_observables_gpu_fields_density
)

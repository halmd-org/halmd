halmd_add_library(halmd_mdsim_gpu_forces
  tabulated_external.cpp
  tabulated_external_kernel.cu
)

halmd_add_modules(
  libhalmd_mdsim_gpu_forces_tabulated_external
)

halmd_add_library(halmd_mdsim_host_integrators
  verlet.cpp
  verlet_nvt_andersen.cpp
)
halmd_add_modules(
  libhalmd_mdsim_host_integrators_verlet
  libhalmd_mdsim_host_integrators_verlet_nvt_andersen
)

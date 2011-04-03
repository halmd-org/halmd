add_library(halmd_mdsim_gpu_integrators_verlet
  verlet_kernel.cu
  verlet.cpp
)
add_library(halmd_mdsim_gpu_integrators_verlet_nvt_andersen
  verlet_nvt_andersen_kernel.cu
  verlet_nvt_andersen.cpp
)

add_library(halmd_mdsim_gpu_forces
  lennard_jones_kernel.cu
  lennard_jones.cpp
  morse_kernel.cu
  morse.cpp
)

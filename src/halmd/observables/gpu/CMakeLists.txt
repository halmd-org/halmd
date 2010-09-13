add_library(halmd_observables_gpu
  thermodynamics.cpp
  thermodynamics_kernel.cu
  )
target_link_libraries(halmd_observables_gpu
  ${CUDA_LIBRARIES}
  )

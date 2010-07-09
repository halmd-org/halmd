add_library(halmd_mdsim_gpu_sort
  hilbert_kernel.cu
  )
target_link_libraries(halmd_mdsim_gpu_sort
  ${CUDA_LIBRARIES}
  )

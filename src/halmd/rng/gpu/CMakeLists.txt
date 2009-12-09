add_library(gpu_rand48
  rand48.cu
  )
target_link_libraries(gpu_rand48
  gpu_prefix_sum
  )

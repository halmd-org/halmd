add_library(ljfluid_base_gpu
  ljfluid_base.cu
  )
add_library(ljfluid_simple_gpu
  ljfluid_simple.cu
  )
add_library(ljfluid_cell_gpu
  ljfluid_cell.cu
  )
add_library(ljfluid_nbr_gpu
  ljfluid_nbr.cu
  )
add_library(hilbert_gpu
  hilbert.cu
  )
add_library(radix_gpu
  radix.cu
  )
add_library(scan_gpu
  scan.cu
  )
add_library(rand48_gpu
  rand48.cu
  scan.cu
  )
set_target_properties(rand48_gpu PROPERTIES
  COMPILE_DEFINITIONS "USE_SCAN_RAND48"
  )

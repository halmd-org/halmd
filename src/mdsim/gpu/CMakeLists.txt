ADD_LIBRARY(ljfluid_base_gpu
  ljfluid_base.cu
  )
ADD_LIBRARY(ljfluid_simple_gpu
  ljfluid_simple.cu
  )
ADD_LIBRARY(ljfluid_cell_gpu
  ljfluid_cell.cu
  )
ADD_LIBRARY(ljfluid_nbr_gpu
  ljfluid_nbr.cu
  )
ADD_LIBRARY(hilbert_gpu
  hilbert.cu
  )
ADD_LIBRARY(radix_gpu
  radix.cu
  )
ADD_LIBRARY(scan_gpu
  scan.cu
  )
ADD_LIBRARY(rand48_gpu
  rand48.cu
  scan.cu
  )
SET_TARGET_PROPERTIES(rand48_gpu PROPERTIES
  COMPILE_DEFINITIONS "USE_SCAN_RAND48"
  )

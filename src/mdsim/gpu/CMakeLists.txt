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

FOREACH(dim RANGE 2 3)
  ADD_LIBRARY(ljfluid${dim}d_base_gpu
    ljfluid_base.cu
    )
  SET_TARGET_PROPERTIES(ljfluid${dim}d_base_gpu PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  ADD_LIBRARY(ljfluid${dim}d_simple_gpu
    ljfluid_simple.cu
    )
  SET_TARGET_PROPERTIES(ljfluid${dim}d_simple_gpu PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  ADD_LIBRARY(ljfluid${dim}d_cell_gpu
    ljfluid_cell.cu
    )
  SET_TARGET_PROPERTIES(ljfluid${dim}d_cell_gpu PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  ADD_LIBRARY(ljfluid${dim}d_nbr_gpu
    ljfluid_nbr.cu
    )
  SET_TARGET_PROPERTIES(ljfluid${dim}d_nbr_gpu PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  ADD_LIBRARY(hilbert${dim}d_gpu
    hilbert.cu
    )
  SET_TARGET_PROPERTIES(hilbert${dim}d_gpu PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
ENDFOREACH(dim)

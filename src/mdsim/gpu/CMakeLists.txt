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
  IF(ljfluid_VARIANT_NEIGHBOUR OR NOT ljfluid_VARIANT_CELL)
    ADD_LIBRARY(ljfluid${dim}d_gpu
      ljfluid_nbr.cu
      )
  ELSE(ljfluid_VARIANT_NEIGHBOUR OR NOT ljfluid_VARIANT_CELL)
    ADD_LIBRARY(ljfluid${dim}d_gpu
      ljfluid_cell.cu
      )
  ENDIF(ljfluid_VARIANT_NEIGHBOUR OR NOT ljfluid_VARIANT_CELL)
  SET_TARGET_PROPERTIES(ljfluid${dim}d_gpu PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )

  IF(ljfluid_VARIANT_NEIGHBOUR OR NOT ljfluid_VARIANT_CELL)
    ADD_LIBRARY(ljfluid${dim}d_smooth_gpu
      ljfluid_nbr.cu
      )
  ELSE(ljfluid_VARIANT_NEIGHBOUR OR NOT ljfluid_VARIANT_CELL)
    ADD_LIBRARY(ljfluid${dim}d_smooth_gpu
      ljfluid_cell.cu
      )
  ENDIF(ljfluid_VARIANT_NEIGHBOUR OR NOT ljfluid_VARIANT_CELL)
  SET_TARGET_PROPERTIES(ljfluid${dim}d_smooth_gpu PROPERTIES
    COMPILE_FLAGS "-DDIM_${dim}D -DUSE_POTENTIAL_SMOOTHING"
    )

  ADD_LIBRARY(hilbert${dim}d_gpu
    hilbert.cu
    )
  SET_TARGET_PROPERTIES(hilbert${dim}d_gpu PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
ENDFOREACH(dim)

ADD_LIBRARY(rand48_gpu rand48.cu)

FOREACH(dim RANGE 2 3)
  ADD_LIBRARY(ljfluid${dim}d_gpu ljfluid.cu)
  SET_TARGET_PROPERTIES(ljfluid${dim}d_gpu PROPERTIES COMPILE_DEFINITIONS "DIM_${dim}D")
ENDFOREACH(dim)

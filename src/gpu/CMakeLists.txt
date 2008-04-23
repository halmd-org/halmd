ADD_LIBRARY(rand48_gpu rand48.cu)

FOREACH(dim RANGE 2 3)
  ADD_LIBRARY(ljfluid${dim}d_gpu ljfluid.cu)
  SET_TARGET_PROPERTIES(ljfluid${dim}d_gpu PROPERTIES COMPILE_FLAGS "-DDIM_${dim}D")

  ADD_LIBRARY(ljfluid${dim}d_leapfrog_gpu ljfluid.cu)
  SET_TARGET_PROPERTIES(ljfluid${dim}d_leapfrog_gpu PROPERTIES COMPILE_FLAGS "-DDIM_${dim}D -DUSE_LEAPFROG")
ENDFOREACH(dim)

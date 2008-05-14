ADD_EXECUTABLE(test_sum sum.cpp sum_gpu.cu)
TARGET_LINK_LIBRARIES(test_sum ${CUDA_LIBRARIES})

FOREACH(dim RANGE 2 3)
  ADD_EXECUTABLE(t_trajectory${dim}d trajectory.cpp)
  SET_TARGET_PROPERTIES(t_trajectory${dim}d PROPERTIES COMPILE_DEFINITIONS "DIM_${dim}D")
  TARGET_LINK_LIBRARIES(t_trajectory${dim}d ${HDF5_CPP_LIBRARY} ${CUDA_LIBRARIES})
ENDFOREACH(dim)

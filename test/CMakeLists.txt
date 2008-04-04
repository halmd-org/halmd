ADD_EXECUTABLE(test_sum sum.cpp sum_gpu.cu)
TARGET_LINK_LIBRARIES(test_sum ${CUDA_LIBRARIES})

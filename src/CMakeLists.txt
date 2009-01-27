add_subdirectory(ljgpu)

if(WITH_CUDA)
  add_subdirectory(cufree)
  add_subdirectory(nvlock)
endif(WITH_CUDA)

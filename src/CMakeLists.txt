add_subdirectory(ljgpu)

if(WITH_CUDA)
  add_subdirectory(cumem)
  add_subdirectory(nvlock)
endif(WITH_CUDA)

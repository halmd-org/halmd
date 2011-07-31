add_subdirectory(host)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

add_subdirectory(halmd)

if(WITH_CUDA AND NOT CUDA_EMULATION)
  add_subdirectory(cufree)
  add_subdirectory(nvlock)
endif(WITH_CUDA AND NOT CUDA_EMULATION)

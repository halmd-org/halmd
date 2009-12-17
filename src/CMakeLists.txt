add_subdirectory(halmd)

if(WITH_CUDA AND NOT CUDA_EMULATION)
  add_subdirectory(cufree)
endif(WITH_CUDA AND NOT CUDA_EMULATION)

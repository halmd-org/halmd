add_executable(cufree cufree.c)
target_link_libraries(cufree
  ${CUDA_LIBRARY}
  )
install(TARGETS cufree
  RUNTIME DESTINATION bin
  )

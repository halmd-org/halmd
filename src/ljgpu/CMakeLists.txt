if(ljgpu_VARIANT_CUDA)
  add_subdirectory(gpu)
endif(ljgpu_VARIANT_CUDA)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Generate version header file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

foreach(dim RANGE 2 3)
  add_executable(ljgpu${dim}d
    H5param.cpp
    H5xx.cpp
    log.cpp
    main.cpp
    options.cpp
    perf.cpp
    )
  target_link_libraries(ljgpu${dim}d
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${HDF5_CPP_LIBRARY}
    )
  if(ljgpu_VARIANT_CUDA)
    target_link_libraries(ljgpu${dim}d
      rand48_gpu
      scan_gpu
      radix_gpu
      ${CUDA_LIBRARIES}
      )
    if(ljgpu_VARIANT_NEIGHBOUR)
      target_link_libraries(ljgpu${dim}d
	ljfluid_nbr_gpu
	hilbert_gpu
	)
    elseif(ljgpu_VARIANT_CELL)
      target_link_libraries(ljgpu${dim}d
	ljfluid_cell_gpu
	)
    else(ljgpu_VARIANT_NEIGHBOUR)
      target_link_libraries(ljgpu${dim}d
	ljfluid_simple_gpu
	)
    endif(ljgpu_VARIANT_NEIGHBOUR)
  else(ljgpu_VARIANT_CUDA)
    target_link_libraries(ljgpu${dim}d
      ${GSL_LIBRARIES}
      rt
      )
  endif(ljgpu_VARIANT_CUDA)
  set_target_properties(ljgpu${dim}d PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  install(TARGETS ljgpu${dim}d
    RUNTIME DESTINATION bin
    )
endforeach(dim)

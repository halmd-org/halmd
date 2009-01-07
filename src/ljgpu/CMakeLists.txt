add_subdirectory(algorithm)
add_subdirectory(mdsim)
add_subdirectory(rng)
add_subdirectory(sample)
add_subdirectory(util)

# Generate version header file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

foreach(dim RANGE 2 3)
  add_executable(ljgpu${dim}d
    main.cpp
    options.cpp
    )
  target_link_libraries(ljgpu${dim}d
    mdsample
    mdutil
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${HDF5_CPP_LIBRARY}
    )
  if(ljgpu_VARIANT_CUDA)
    target_link_libraries(ljgpu${dim}d
      gpu_lattice
      gpu_rand48
      ${CUDA_LIBRARIES}
      )
    if(ljgpu_VARIANT_NEIGHBOUR)
      target_link_libraries(ljgpu${dim}d
	gpu_hilbert
	gpu_ljfluid_nbr
	gpu_radix_sort
	)
    elseif(ljgpu_VARIANT_CELL)
      target_link_libraries(ljgpu${dim}d
	gpu_ljfluid_cell
	)
    else(ljgpu_VARIANT_NEIGHBOUR)
      target_link_libraries(ljgpu${dim}d
	gpu_ljfluid_square
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

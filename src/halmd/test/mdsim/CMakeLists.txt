if(WITH_CUDA)
  add_executable(sample_smooth_function smooth.cpp)
  set_target_properties(sample_smooth_function PROPERTIES
    COMPILE_FLAGS "-DDIM_3D"
    )
  target_link_libraries(sample_smooth_function
    gpu_ljfluid_square
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    pthread
    )

  foreach(dim RANGE 2 3)
    add_executable(hilbert_curve_sort${dim}d hilbert.cpp)
    set_target_properties(hilbert_curve_sort${dim}d PROPERTIES
      COMPILE_FLAGS "-DDIM_${dim}D"
      )
    target_link_libraries(hilbert_curve_sort${dim}d
      gpu_hilbert
      gpu_lattice
      gpu_radix_sort
      gpu_rand48
      ${Boost_PROGRAM_OPTIONS_LIBRARY}
      ${CUDA_LIBRARIES}
      pthread
      rt
      )
  endforeach(dim)
endif(WITH_CUDA)

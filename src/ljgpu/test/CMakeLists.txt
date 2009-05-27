if(WITH_CUDA)
  add_subdirectory(dsfun)

  add_executable(sample_smooth_function smooth.cpp)
  set_target_properties(sample_smooth_function PROPERTIES
    COMPILE_FLAGS "-DDIM_3D"
    )
  target_link_libraries(sample_smooth_function
    gpu_ljfluid_square
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )

  add_executable(parallel_scan scan.cpp)
  target_link_libraries(parallel_scan
    gpu_prefix_sum
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )
  foreach(count 10000000 1000000 100000 10000 1024 1000 512 100 10 2)
    add_test(parallel_scan_${count} parallel_scan -N${count})
  endforeach(count)

  add_executable(parallel_radix_sort radix.cpp)
  target_link_libraries(parallel_radix_sort
    gpu_radix_sort
    gpu_rand48
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
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
      )
  endforeach(dim)

  add_executable(rand48_test rand48.cpp)
  target_link_libraries(rand48_test
    gpu_rand48
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )
endif(WITH_CUDA)

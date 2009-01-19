ENABLE_TESTING()

ADD_EXECUTABLE(sample_smooth_function smooth.cpp)
SET_TARGET_PROPERTIES(sample_smooth_function PROPERTIES
  COMPILE_FLAGS "-DDIM_3D"
  )
TARGET_LINK_LIBRARIES(sample_smooth_function
  gpu_ljfluid_square
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${CUDA_LIBRARIES}
  )

ADD_EXECUTABLE(parallel_scan scan.cpp)
TARGET_LINK_LIBRARIES(parallel_scan
  gpu_prefix_sum
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${CUDA_LIBRARIES}
  )
FOREACH(count 10000000 1000000 100000 10000 1024 1000 512 100 10 2)
  ADD_TEST(parallel_scan_${count} parallel_scan -N${count})
ENDFOREACH(count)

ADD_EXECUTABLE(parallel_radix_sort radix.cpp)
TARGET_LINK_LIBRARIES(parallel_radix_sort
  gpu_radix_sort
  gpu_rand48
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${CUDA_LIBRARIES}
  )

FOREACH(dim RANGE 2 3)
  ADD_EXECUTABLE(hilbert_curve_sort${dim}d hilbert.cpp)
  SET_TARGET_PROPERTIES(hilbert_curve_sort${dim}d PROPERTIES
    COMPILE_FLAGS "-DDIM_${dim}D"
    )
  TARGET_LINK_LIBRARIES(hilbert_curve_sort${dim}d
    gpu_hilbert
    gpu_lattice
    gpu_radix_sort
    gpu_rand48
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )
ENDFOREACH(dim)

ADD_EXECUTABLE(rand48_test rand48.cpp)
TARGET_LINK_LIBRARIES(rand48_test
  gpu_rand48
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${CUDA_LIBRARIES}
  )

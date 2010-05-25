if(WITH_CUDA)
  add_executable(sample_smooth_function smooth.cpp)
  set_target_properties(sample_smooth_function PROPERTIES
    COMPILE_FLAGS "-DDIM_3D"
    )
  target_link_libraries(sample_smooth_function
    gpu_ljfluid_square
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
endif(WITH_CUDA)

foreach(impl ${HALMD_impl})
  string(REGEX REPLACE "_impl[_a-z]*$" "" class ${impl})
  string(REGEX REPLACE "(^ljfluid_impl_|_impl)" "" backend ${impl})

  add_executable(test_thermodynamics_${backend}
    thermodynamics.cpp
    ../../options.cpp
    ../../mdlib.cpp
  )
  set(HALMD_COMPILE_DEFINITIONS
    MDSIM_CLASS=${class}
    MDSIM_IMPL=${impl}
    MDSIM_BACKEND="${backend}"
    )
  add_definitions(-DBACKEND_EXECUTABLES)

  target_link_libraries(test_thermodynamics_${backend}
    ${HALMD_COMMON_LIBRARIES}
  )
  if(${backend} STREQUAL "host")
    target_link_libraries(test_thermodynamics_${backend}
      ${GSL_LIBRARIES}
      sample_host
    )
    set_target_properties(test_thermodynamics_${backend} PROPERTIES
      COMPILE_FLAGS "-UWITH_CUDA"
    )
  else(${backend} STREQUAL "host")
    target_link_libraries(test_thermodynamics_${backend}
      ${CUDA_LIBRARIES}
      gpu_boltzmann
      gpu_lattice
      gpu_radix_sort
      gpu_rand48
      gpu_reduce
      gpu_virial
      sample_gpu
    )
    if(${backend} STREQUAL "gpu_neighbour")
      target_link_libraries(test_thermodynamics_${backend}
        gpu_ljfluid_nbr
        gpu_hilbert
      )
    endif(${backend} STREQUAL "gpu_neighbour")
    if(${backend} STREQUAL "gpu_square")
      target_link_libraries(test_thermodynamics_${backend}
        gpu_ljfluid_square
      )
    endif(${backend} STREQUAL "gpu_square")
  endif(${backend} STREQUAL "host")
  set_target_properties(test_thermodynamics_${backend} PROPERTIES
    COMPILE_DEFINITIONS "${HALMD_COMPILE_DEFINITIONS}"
    )

  add_test(thermodynamics_${backend} test_thermodynamics_${backend})
endforeach(impl)


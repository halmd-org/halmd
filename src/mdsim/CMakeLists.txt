IF(ljfluid_VARIANT_CUDA)
  ADD_SUBDIRECTORY(gpu)
ENDIF(ljfluid_VARIANT_CUDA)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Generate version header file
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

FOREACH(dim RANGE 2 3)
  ADD_EXECUTABLE(ljfluid${dim}d
    H5param.cpp
    H5xx.cpp
    log.cpp
    main.cpp
    options.cpp
    perf.cpp
    )
  TARGET_LINK_LIBRARIES(ljfluid${dim}d
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${HDF5_CPP_LIBRARY}
    )
  IF(ljfluid_VARIANT_CUDA)
    TARGET_LINK_LIBRARIES(ljfluid${dim}d
      rand48_gpu
      scan_gpu
      radix_gpu
      ${CUDA_LIBRARIES}
      )
    IF(ljfluid_VARIANT_NEIGHBOUR)
      TARGET_LINK_LIBRARIES(ljfluid${dim}d
	ljfluid${dim}d_nbr_gpu
	hilbert${dim}d_gpu
	)
    ELSEIF(ljfluid_VARIANT_CELL)
      TARGET_LINK_LIBRARIES(ljfluid${dim}d
	ljfluid${dim}d_cell_gpu
	)
    ELSE(ljfluid_VARIANT_NEIGHBOUR)
      TARGET_LINK_LIBRARIES(ljfluid${dim}d
	ljfluid${dim}d_simple_gpu
	)
    ENDIF(ljfluid_VARIANT_NEIGHBOUR)
  ELSE(ljfluid_VARIANT_CUDA)
    TARGET_LINK_LIBRARIES(ljfluid${dim}d
      ${GSL_LIBRARIES}
      rt
      )
  ENDIF(ljfluid_VARIANT_CUDA)
  SET_TARGET_PROPERTIES(ljfluid${dim}d PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  INSTALL(TARGETS ljfluid${dim}d
    RUNTIME DESTINATION bin
    )
ENDFOREACH(dim)

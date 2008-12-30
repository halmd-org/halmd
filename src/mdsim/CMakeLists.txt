# Add program variant compile flags to cmake GUI
SET(ljfluid_VARIANT_CUDA TRUE CACHE BOOL
  "Use CUDA")
IF(ljfluid_VARIANT_CUDA)
  ADD_DEFINITIONS(-DUSE_CUDA)
ENDIF(ljfluid_VARIANT_CUDA)

SET(ljfluid_VARIANT_CELL CACHE BOOL
  "Use cell lists and Verlet neighbour lists")
IF(ljfluid_VARIANT_CELL)
  ADD_DEFINITIONS(-DUSE_CELL)
ENDIF(ljfluid_VARIANT_CELL)

SET(ljfluid_VARIANT_HILBERT CACHE BOOL
  "Use Hilbert space-filling curve particle ordering")
IF(ljfluid_VARIANT_HILBERT)
  ADD_DEFINITIONS(-DUSE_HILBERT_ORDER)
ENDIF(ljfluid_VARIANT_HILBERT)

SET(ljfluid_VARIANT_C2POT CACHE BOOL
  "Use C2 potential smoothing")
IF(ljfluid_VARIANT_DSFUN)
  ADD_DEFINITIONS(-DUSE_DSFUN)
ENDIF(ljfluid_VARIANT_DSFUN)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Generate version header file
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

FOREACH(dim RANGE 2 3)
  IF(ljfluid_VARIANT_CUDA)
    ADD_EXECUTABLE(ljfluid${dim}d
      H5param.cpp
      H5xx.cpp
      correlation.cpp
      energy.cpp
      ljfluid_gpu.cpp
      log.cpp
      main.cpp
      mdsim.cpp
      options.cpp
      perf.cpp
      trajectory.cpp
      )
  ELSE(ljfluid_VARIANT_CUDA)
    ADD_EXECUTABLE(ljfluid${dim}d
      H5param.cpp
      H5xx.cpp
      correlation.cpp
      energy.cpp
      ljfluid_host.cpp
      log.cpp
      main.cpp
      mdsim.cpp
      options.cpp
      perf.cpp
      trajectory.cpp
      )
  ENDIF(ljfluid_VARIANT_CUDA)
  SET_TARGET_PROPERTIES(ljfluid${dim}d PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  IF(ljfluid_VARIANT_CUDA)
    TARGET_LINK_LIBRARIES(ljfluid${dim}d
      ljfluid${dim}d_gpu
      hilbert${dim}d_gpu
      rand48_gpu
      scan_gpu
      radix_gpu
      ${Boost_FILESYSTEM_LIBRARY}
      ${Boost_PROGRAM_OPTIONS_LIBRARY}
      ${CUDA_LIBRARIES}
      ${HDF5_CPP_LIBRARY}
      rt
      )
  ELSE(ljfluid_VARIANT_CUDA)
    TARGET_LINK_LIBRARIES(ljfluid${dim}d
      ${Boost_FILESYSTEM_LIBRARY}
      ${Boost_PROGRAM_OPTIONS_LIBRARY}
      ${GSL_LIBRARIES}
      ${HDF5_CPP_LIBRARY}
      rt
      )
  ENDIF(ljfluid_VARIANT_CUDA)
  INSTALL(TARGETS ljfluid${dim}d
    RUNTIME DESTINATION bin
    )
ENDFOREACH(dim)

IF(ljfluid_VARIANT_CUDA)
  ADD_SUBDIRECTORY(gpu)
ENDIF(ljfluid_VARIANT_CUDA)

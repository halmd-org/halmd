# Add program variant compile flags to cmake GUI
SET(ljfluid_VARIANT_C2POT CACHE BOOL
  "Use C2 potential smoothing")
IF(ljfluid_VARIANT_DSFUN)
  ADD_DEFINITIONS(-DUSE_DSFUN)
ENDIF(ljfluid_VARIANT_DSFUN)

SET(ljfluid_VARIANT_CELL CACHE BOOL
  "Use cell lists")
IF(ljfluid_VARIANT_CELL)
  ADD_DEFINITIONS(-DUSE_CELL)
ENDIF(ljfluid_VARIANT_CELL)

SET(ljfluid_VARIANT_CELL_SUMMATION_ORDER CACHE BOOL
  "Use opposite cell summation order")
IF(ljfluid_VARIANT_CELL_SUMMATION_ORDER)
  ADD_DEFINITIONS(-DUSE_CELL_SUMMATION_ORDER)
ENDIF(ljfluid_VARIANT_CELL_SUMMATION_ORDER)

SET(ljfluid_VARIANT_DSFUN CACHE BOOL
  "Use double-single precision functions")
IF(ljfluid_VARIANT_C2POT)
  ADD_DEFINITIONS(-DUSE_POTENTIAL_SMOOTHING)
ENDIF(ljfluid_VARIANT_C2POT)

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
    correlation.cpp
    energy.cpp
    ljfluid.cpp
    log.cpp
    main.cpp
    mdsim.cpp
    options.cpp
    perf.cpp
    trajectory.cpp
    )
  SET_TARGET_PROPERTIES(ljfluid${dim}d PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  TARGET_LINK_LIBRARIES(ljfluid${dim}d
    ljfluid${dim}d_gpu
    rand48_gpu
    scan_gpu
    radix_gpu
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${GSL_LIBRARIES}
    ${CUDA_LIBRARIES}
    ${HDF5_CPP_LIBRARY}
    rt
    )
  INSTALL(TARGETS ljfluid${dim}d
    RUNTIME DESTINATION bin
    )
ENDFOREACH(dim)

ADD_SUBDIRECTORY(gpu)

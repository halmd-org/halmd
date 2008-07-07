INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Generate version header file
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

FOREACH(dim RANGE 2 3)
  ADD_EXECUTABLE(ljfluid${dim}d main.cpp options.cpp H5param.cpp H5xx.cpp log.cpp)
  SET_TARGET_PROPERTIES(ljfluid${dim}d PROPERTIES COMPILE_DEFINITIONS "DIM_${dim}D")
  TARGET_LINK_LIBRARIES(ljfluid${dim}d
    ljfluid${dim}d_gpu
    rand48_gpu
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${GSL_LIBRARIES}
    ${CUDA_LIBRARIES}
    ${HDF5_CPP_LIBRARY}
    )
ENDFOREACH(dim)

ADD_SUBDIRECTORY(gpu)


ADD_DEFINITIONS(-Wall)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Generate version header file
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )


FOREACH(dim RANGE 2 3)
  ADD_EXECUTABLE(hardspheres${dim}d main.cpp options.cpp H5param.cpp H5xx.cpp log.cpp)
  SET_TARGET_PROPERTIES(hardspheres${dim}d PROPERTIES COMPILE_DEFINITIONS "DIM_${dim}D")
  TARGET_LINK_LIBRARIES(hardspheres${dim}d
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${GSL_LIBRARIES}
    ${HDF5_CPP_LIBRARY}
    rt
    )
ENDFOREACH(dim)

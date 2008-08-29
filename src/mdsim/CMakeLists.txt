ADD_DEFINITIONS(-Wall)

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
    options.cpp
    perf.cpp
    )
  SET_TARGET_PROPERTIES(ljfluid${dim}d PROPERTIES
    COMPILE_DEFINITIONS "DIM_${dim}D"
    )
  TARGET_LINK_LIBRARIES(ljfluid${dim}d
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${GSL_LIBRARIES}
    ${HDF5_CPP_LIBRARY}
    rt
    )
ENDFOREACH(dim)

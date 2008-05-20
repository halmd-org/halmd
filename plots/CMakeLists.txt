FOREACH(filename ljplot trajectories2d.gpi trajectories3d.gpi gnuplot_movie)
  CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/${filename}"
    "${CMAKE_CURRENT_BINARY_DIR}/${filename}"
    @ONLY)
ENDFOREACH(filename)

FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

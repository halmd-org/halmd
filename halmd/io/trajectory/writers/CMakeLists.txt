add_library(halmd_io_trajectory_writers
  hdf5.cpp
  )
target_link_libraries(halmd_io_trajectory_writers
  ${HDF5_CPP_LIBRARY}
  ${HDF5_LIBRARY}
  )

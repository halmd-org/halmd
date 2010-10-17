add_library(halmd_io_trajectory_readers
  hdf5.cpp
  )
target_link_libraries(halmd_io_trajectory_readers
  ${HDF5_CPP_LIBRARY}
  ${HDF5_LIBRARY}
  )

if(CMAKE_CXX_PLATFORM_ID STREQUAL "Linux")

  # On Linux, add --Wl,--as-needed to the default linker flags.
  # This makes distribution packagers happy, as binaries will only
  # reference the actually used libraries, which allows for minimal package
  # dependencies of the halmd package. The HALMD_COMMON_LIBRARIES variable
  # contains more libraries than needed in most cases, to support linking
  # Boost and HDF5 libraries statically.

  set(CMAKE_EXE_LINKER_FLAGS_INIT "-Wl,--as-needed")
  set(CMAKE_MODULE_LINKER_FLAGS_INIT "-Wl,--as-needed")
  set(CMAKE_SHARED_LINKER_FLAGS_INIT "-Wl,--as-needed")
endif()

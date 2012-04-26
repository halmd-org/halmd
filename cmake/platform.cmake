if(DEFINED CMAKE_CXX_COMPILER_ID)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall -std=c++98")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG -fvisibility-inlines-hidden")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "XL")

    set(CMAKE_CXX_FLAGS_INIT "-qrtti=all")

    # We only need object-level inlining (opposed to link-time optimization),
    # therefore we follow the recommendation for XL 11 that -qipa=inline is
    # deprecated, and use -qinline instead.
    # http://publib.boulder.ibm.com/infocenter/iadthelp/v8r0/index.jsp?topic=/com.ibm.xlcpp111.linux.doc/compiler_ref/opt_ipa.html
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O4 -DNDEBUG -qstrict=all -qnoipa -qinline")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O4 -DNDEBUG -qstrict=all -qnoipa -qinline")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-g -O2 -qstrict=all")

  else()
    message(WARNING "Unsupported CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
  endif()
endif()

if(DEFINED CMAKE_CUDA_COMPILER_ID)
  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVCC")

    set(CMAKE_CUDA_FLAGS_INIT "-Xcompiler -fPIC -Xptxas -v -arch sm_12")

  else()
    message(WARNING "Unsupported CUDA compiler: ${CMAKE_CUDA_COMPILER_ID}")
  endif()
endif()

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

  # Strip binaries for Release builds
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE_INIT "-Wl,-s")
  set(CMAKE_MODULE_LINKER_FLAGS_RELEASE_INIT "-Wl,-s")
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE_INIT "-Wl,-s")

elseif(CMAKE_CXX_PLATFORM_ID STREQUAL "AIX")

  # The combined flags -qtwolink -Wl,-bweaklocal instruct the AIX linker
  # to discard unreferenced objects of an archive, which is needed to
  # avoid undefined symbol errors when linking the unit tests.
  # -qnoipa is needed to disable link-time IPA, which breaks -qtwolink.
  set(CMAKE_EXE_LINKER_FLAGS_INIT "-qtwolink -Wl,-bweaklocal -qnoipa")
  set(CMAKE_MODULE_LINKER_FLAGS_INIT "-qtwolink -Wl,-bweaklocal -qnoipa")
  set(CMAKE_SHARED_LINKER_FLAGS_INIT "-qtwolink -Wl,-bweaklocal -qnoipa")
endif()

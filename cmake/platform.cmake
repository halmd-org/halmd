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

    # Compile for CUDA compute version 1.2, and generate PTX 1.2 code
    # as well as binary code for targets of compute capability 1.3 and 2.0.
    # This ensures backward and forward compatibility with targets of
    # compute capability ≥ 1.2 through JIT compilation at program
    # startup, while providing binary code for Tesla and Fermi GPUs.
    #
    # Note that we compile with -arch=compute_12 to disable native double
    # precision (present with compute_13), and to avoid performance penalty
    # of IEEE-compliant floating-point (present with compute_20).
    #
    # We will raise the default compute version later to enable native double
    # precision, as an alternative to double-single precision, and when all
    # kernels are optimised for Fermi GPUs.

    set(CMAKE_CUDA_FLAGS_INIT "-Xcompiler -fPIC -Xptxas -v -arch=compute_12 -code=compute_12,sm_13,sm_20")

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

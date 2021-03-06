set(CMAKE_BUILD_TYPE_INIT "Release")

if(DEFINED CMAKE_CXX_COMPILER_ID)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    # HALMD requires a C++11 compiler, e.g. GCC 4.7.
    # Remove -DNDEBUG from RelWithDebInfo to enable assert() and LOG_DEBUG/LOG_TRACE.
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
      message(FATAL_ERROR "Minimal supported version of GCC compiler is 4.7")
    endif()
    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall -std=c++11 -pedantic")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS -fvisibility=hidden")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

    # Clang versions >= 3.9.1 issue warnings if a template specialization is
    # declared, but not defined in a compilation unit.
    # HALMD regularly uses explicit template specializations in specific compilation
    # units. This paradigm conforms to the C++ standards and is intended,
    # therefore the warnings are disabled.
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.9.1")
      set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall -std=c++11 -pedantic -Wno-undefined-var-template")
    else()
      set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall -std=c++11 -pedantic")
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS -fvisibility=hidden")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
    # clang doesn't print colored diagnostics when invoked from Ninja
    if (UNIX AND CMAKE_GENERATOR STREQUAL "Ninja")
      set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -fcolor-diagnostics")
    endif()

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

    set(CMAKE_CXX_FLAGS_INIT "-fPIC -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS -fvisibility=hidden")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

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

    # Compile for CUDA compute capability 3.5 (Kepler), and generate PTX 3.5 code
    # as well as binary code for this target. (The sm_XX parameter is adjusted
    # best in the build tree after the initial configuration.) This ensures
    # backward and forward compatibility with targets of compute capability ≥
    # 1.2 through JIT compilation at program startup, while providing binary
    # code for Kepler GPUs.
    #
    # Note that the MD core modules were developed for compute capability 1.2,
    # which brings only non-IEEE-compliant, single-precision floating-point
    # arithmetics. So we consider it safe to disable IEEE-compliance, which
    # has (small) performance penalties.
    #
    # CUDA versions starting from 6.5 do not support comute version 1.2 anymore,
    # and 2.0 was deprecated in CUDA 9.0, so in these cases the default compute
    # version is raised to 3.5 (Kepler) regardless.

    #
    # the following flags must be duplicated in /CMakeLists.txt for CMake without native CUDA support
    #
    if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 6.5)
      set(CMAKE_CUDA_FLAGS_INIT "-Xcompiler -fPIC -Xptxas -v -std=c++11 -arch=compute_12 -code=compute_12")
    else()
      set(CMAKE_CUDA_FLAGS_INIT "-Xcompiler -fPIC -Xptxas -v -std=c++11 -arch=compute_35 -code=compute_35,sm_35 -ftz=true -prec-div=false -prec-sqrt=false --fmad=true")
    endif()

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

  if(HALMD_USE_STATIC_LIBS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.5")
      if(NOT DEFINED HALMD_USE_STATIC_LIBSTDCXX)
        set(HALMD_USE_STATIC_LIBSTDCXX TRUE)
      endif()
      if(NOT DEFINED HALMD_USE_STATIC_LIBGCC)
        set(HALMD_USE_STATIC_LIBGCC TRUE)
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.1")
      if(NOT DEFINED HALMD_USE_STATIC_LIBSTDCXX)
        set(HALMD_USE_STATIC_LIBSTDCXX TRUE)
      endif()
      if(NOT DEFINED HALMD_USE_STATIC_LIBGCC)
        set(HALMD_USE_STATIC_LIBGCC TRUE)
      endif()
    endif()
  endif()

  if(HALMD_USE_STATIC_LIBSTDCXX)
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -static-libstdc++")
  endif()
  if(HALMD_USE_STATIC_LIBGCC)
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -static-libgcc")
  endif()

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

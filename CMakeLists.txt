cmake_minimum_required(VERSION 2.8)

if(NOT HALMD_DOC_ONLY)
  set(HALMD_DOC_ONLY FALSE)
else()
  set(HALMD_DOC_ONLY TRUE)
endif()
set(HALMD_DOC_ONLY ${HALMD_DOC_ONLY} CACHE INTERNAL "Build HALMD documentation only" FORCE)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Git QUIET REQUIRED)
find_package(GitRepository QUIET REQUIRED)

git_repository("${CMAKE_SOURCE_DIR}" HALMD)

set(PROGRAM_NAME "halmd")
set(PROGRAM_DESC "HAL’s MD package")
set(PROGRAM_AUTHORS "Peter Colberg and Felix Höfling")
set(PROGRAM_COPYRIGHT "2007-2012, ${PROGRAM_AUTHORS}")
if(HALMD_GIT_COMMIT_TAG)
  string(REGEX MATCH "v[0-9].*" PROGRAM_VERSION "${HALMD_GIT_COMMIT_TAG}")
  if(NOT PROGRAM_VERSION)
    set(PROGRAM_VERSION "${HALMD_GIT_COMMIT_TAG}")
    if(HALMD_GIT_BRANCH)
      set(PROGRAM_VERSION "${PROGRAM_VERSION} (${HALMD_GIT_BRANCH})")
    endif()
  endif()
else()
  set(PROGRAM_VERSION "(unknown version)")
endif()
if(HALMD_GIT_COMMITTER_DATE)
  set(PROGRAM_DATE "${HALMD_GIT_COMMITTER_DATE}")
else()
  set(PROGRAM_DATE "unknown date")
endif()

# Output source version for Dashboard test reports
message(STATUS "Building HALMD ${PROGRAM_VERSION}")

if(NOT HALMD_DOC_ONLY)

project(HALMD CXX)

if(CMAKE_BUILD_TYPE STREQUAL "${CMAKE_BUILD_TYPE_INIT}")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if (HALMD_USE_STATIC_LIBS)
  set(HALMD_USE_STATIC_LIBS TRUE CACHE INTERNAL
    "Build statically linked backend executables")
  set(Boost_USE_STATIC_LIBS TRUE)
  set(HDF5_USE_STATIC_LIBS TRUE)
endif(HALMD_USE_STATIC_LIBS)

# The FindBoost CMake module prefers multi-threaded libraries (filenames with
# postfix "-mt") over non-multi-threaded libraries. On Redhat or SuSE with
# installed system Boost libraries, this causes the system libraries (with
# "-mt") to override the custom-compiled libraries (without "-mt"). For now,
# HALMD requires custom-compiled libraries in any case, as upstream Boost does
# not contain the Boost.Log library yet.
#
if(NOT DEFINED Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED FALSE)
endif(NOT DEFINED Boost_USE_MULTITHREADED)

# Use Luabind debug library if build type is not Release.
# The debug variant has many assert() statements, e.g. to check
# for duplicate class registration, or missing base classes.
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
  set(LUABIND_USE_DEBUG_LIBS TRUE)
endif()

find_package(Boost 1.43.0 QUIET REQUIRED COMPONENTS
  date_time
  filesystem
  log
  program_options
  random
  system
  thread
  unit_test_framework
)
find_package(CUDA QUIET)
find_package(HDF5 QUIET REQUIRED)
find_package(Lua QUIET REQUIRED)
find_package(Luabind QUIET REQUIRED)

if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|Intel)$")
  if(CMAKE_CXX_FLAGS STREQUAL "${CMAKE_CXX_FLAGS_INIT}")
    set(CMAKE_CXX_FLAGS "-fPIC -Wall -std=c++98" CACHE STRING
      "Flags used by the compiler during all build types." FORCE)
  endif()
  if(CMAKE_CXX_FLAGS_RELEASE STREQUAL "${CMAKE_CXX_FLAGS_RELEASE_INIT}")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fvisibility-inlines-hidden" CACHE STRING
      "Flags used by the compiler during release builds." FORCE)
  endif()
  # Strip binaries for Release builds
  if(CMAKE_EXE_LINKER_FLAGS_RELEASE STREQUAL "${CMAKE_EXE_LINKER_FLAGS_RELEASE_INIT}")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-Wl,-s" CACHE STRING
      "Flags used by the linker during release builds." FORCE)
  endif()
  if(CMAKE_MODULE_LINKER_FLAGS_RELEASE STREQUAL "${CMAKE_MODULE_LINKER_FLAGS_RELEASE_INIT}")
    set(CMAKE_MODULE_LINKER_FLAGS_RELEASE "-Wl,-s" CACHE STRING
      "Flags used by the linker during release builds." FORCE)
  endif()
  if(CMAKE_SHARED_LINKER_FLAGS_RELEASE STREQUAL "${CMAKE_SHARED_LINKER_FLAGS_RELEASE_INIT}")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "-Wl,-s" CACHE STRING
      "Flags used by the linker during release builds." FORCE)
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "XL")
  if(CMAKE_CXX_FLAGS STREQUAL "${CMAKE_CXX_FLAGS_INIT}")
    set(CMAKE_CXX_FLAGS "-qrtti=all" CACHE STRING
      "Flags used by the compiler during all build types." FORCE)
  endif()
  # We only need object-level inlining (opposed to link-time optimization),
  # therefore we follow the recommendation for XL 11 that -qipa=inline is
  # deprecated, and use -qinline instead.
  # http://publib.boulder.ibm.com/infocenter/iadthelp/v8r0/index.jsp?topic=/com.ibm.xlcpp111.linux.doc/compiler_ref/opt_ipa.html
  if(CMAKE_CXX_FLAGS_RELEASE STREQUAL "${CMAKE_CXX_FLAGS_RELEASE_INIT}")
    set(CMAKE_CXX_FLAGS_RELEASE "-O4 -DNDEBUG -qstrict=all -qnoipa -qinline" CACHE STRING
      "Flags used by the compiler during Release builds" FORCE)
  endif()
  if(CMAKE_CXX_FLAGS_MINSIZEREL STREQUAL "${CMAKE_CXX_FLAGS_MINSIZEREL_INIT}")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING
      "Flags used by the compiler during MinSizeRel builds" FORCE)
  endif()
  if(CMAKE_CXX_FLAGS_RELWITHDEBINFO STREQUAL "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -qstrict=all" CACHE STRING
      "Flags used by the compiler during RelWithDebInfo builds" FORCE)
  endif()
else()
  message(WARNING "Unsupported CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

if(CMAKE_CXX_PLATFORM_ID STREQUAL "AIX")
  # The combined flags -qtwolink -Wl,-bweaklocal instruct the AIX linker
  # to discard unreferenced objects of an archive, which is needed to
  # avoid undefined symbol errors when linking the unit tests.
  # -qnoipa is needed to disable link-time IPA, which breaks -qtwolink.
  if(CMAKE_EXE_LINKER_FLAGS STREQUAL "${CMAKE_EXE_LINKER_FLAGS_INIT}")
    set(CMAKE_EXE_LINKER_FLAGS "-qtwolink -Wl,-bweaklocal -qnoipa" CACHE STRING
      "Flags used by the linker." FORCE)
  endif()
  if(CMAKE_SHARED_LINKER_FLAGS STREQUAL "${CMAKE_SHARED_LINKER_FLAGS_INIT}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING
      "Flags used by the linker during the creation of shared libraries." FORCE)
  endif()
endif()

if(CUDA_FOUND AND NOT HALMD_USE_STATIC_LIBS)
  set(WITH_CUDA TRUE)
  enable_language(CUDA)

  if(CMAKE_CUDA_FLAGS STREQUAL "${CMAKE_CUDA_FLAGS_INIT}")
    set(CMAKE_CUDA_FLAGS "-Xcompiler -fPIC -Xptxas -v -arch sm_12" CACHE STRING
      "Flags used by the compiler during all build types." FORCE)
  endif()

  add_definitions(-DWITH_CUDA)
else(CUDA_FOUND AND NOT HALMD_USE_STATIC_LIBS)
  set(WITH_CUDA FALSE)
endif(CUDA_FOUND AND NOT HALMD_USE_STATIC_LIBS)

if(HALMD_USE_STATIC_LIBS)
  set(HALMD_BACKEND_EXECUTABLES TRUE CACHE BOOL
    "Build separate executable for each backend" FORCE)
else(HALMD_USE_STATIC_LIBS)
  set(HALMD_BACKEND_EXECUTABLES FALSE CACHE BOOL
    "Build separate executable for each backend")
endif(HALMD_USE_STATIC_LIBS)
if(HALMD_BACKEND_EXECUTABLES)
  add_definitions(-DBACKEND_EXECUTABLES)
endif(HALMD_BACKEND_EXECUTABLES)

set(HALMD_VARIANT_HILBERT_ORDER TRUE CACHE BOOL
  "Use Hilbert space-filling curve particle ordering")
if(HALMD_VARIANT_HILBERT_ORDER)
  add_definitions(-DUSE_HILBERT_ORDER)
endif(HALMD_VARIANT_HILBERT_ORDER)

set(HALMD_VARIANT_HILBERT_ALT_3D FALSE CACHE BOOL
  "Use alternative 3D Hilbert curve vertex rules")
if(HALMD_VARIANT_HILBERT_ALT_3D)
  add_definitions(-DUSE_HILBERT_ALT_3D)
endif(HALMD_VARIANT_HILBERT_ALT_3D)

if(WITH_CUDA)
  set(HALMD_VARIANT_CELL_SUMMATION_ORDER TRUE CACHE BOOL
    "Use opposite cell summation order")
  if(HALMD_VARIANT_CELL_SUMMATION_ORDER)
    add_definitions(-DUSE_CELL_SUMMATION_ORDER)
  endif(HALMD_VARIANT_CELL_SUMMATION_ORDER)

  set(HALMD_VARIANT_FORCE_DSFUN TRUE CACHE BOOL
    "Use double-single precision functions in cell summation")
  if(HALMD_VARIANT_FORCE_DSFUN)
    add_definitions(-DUSE_FORCE_DSFUN)
  endif(HALMD_VARIANT_FORCE_DSFUN)

  set(HALMD_VARIANT_VERLET_DSFUN TRUE CACHE BOOL
    "Use double-single precision functions in Verlet integrator")
  if(HALMD_VARIANT_VERLET_DSFUN)
    add_definitions(-DUSE_VERLET_DSFUN)
  endif(HALMD_VARIANT_VERLET_DSFUN)

  set(HALMD_DEVICE_SCALE "3" CACHE STRING
    "Scale/size of the CUDA device (try to reduce in case of insufficient resources)")
  add_definitions(-DDEVICE_SCALE=${HALMD_DEVICE_SCALE})

endif(WITH_CUDA)

#
# The following option only works on x86-64 by default, which always uses the
# SSE instruction set for floating-point math. On i386, the x87 floating-point
# unit provides 80-bit extended double precision math internally, causing
# excess precision even if values are stored in single-precision.
#
# See the -mfpmath option in the gcc manpage for details, and
#
# Deterministic cross-platform floating point arithmetics
# http://www.christian-seiler.de/projekte/fpmath/
#
set(HALMD_VARIANT_HOST_SINGLE_PRECISION FALSE CACHE BOOL
  "Use single-precision math in host implementation (requires SSE)")
if(HALMD_VARIANT_HOST_SINGLE_PRECISION)
  add_definitions(-DUSE_HOST_SINGLE_PRECISION)
endif(HALMD_VARIANT_HOST_SINGLE_PRECISION)

# set appropriate RPATH on installed binaries as well as in build tree
#
# see http://www.vtk.org/Wiki/CMake_RPATH_handling
#
# use, i.e. don't skip, the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(HALMD_COMMON_LIBRARIES
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_LOG_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_RANDOM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  ${HDF5_CPP_LIBRARY}
  ${HDF5_LIBRARY}
  ${LUABIND_LIBRARY}
  ${LUA_LIBRARIES}
)
if(WITH_CUDA)
  list(APPEND HALMD_COMMON_LIBRARIES
    ${CUDA_LIBRARIES}
  )
endif(WITH_CUDA)
if(HALMD_USE_STATIC_LIBS)
  list(APPEND HALMD_COMMON_LIBRARIES
    -static-libgcc
    -Wl,-Bstatic
  )
endif(HALMD_USE_STATIC_LIBS)
list(APPEND HALMD_COMMON_LIBRARIES
  rt
  dl
  pthread
  z
)

enable_testing()
include(CTest)

include_directories(${HALMD_SOURCE_DIR}/include)
include_directories(${HALMD_SOURCE_DIR}/libs/h5xx)
include_directories(${HALMD_SOURCE_DIR})
include_directories(${HALMD_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${LUA_INCLUDE_DIR})
include_directories(${LUABIND_INCLUDE_DIR})
if(WITH_CUDA)
  include_directories(${CUDA_INCLUDE_DIR})
endif(WITH_CUDA)

add_subdirectory(doc EXCLUDE_FROM_ALL)
add_subdirectory(halmd)
add_subdirectory(lua)
add_subdirectory(test)
add_subdirectory(libs/h5xx/test)

else(NOT HALMD_DOC_ONLY)

project(HALMD NONE) # no compiler languages

add_subdirectory(doc)

endif(NOT HALMD_DOC_ONLY)

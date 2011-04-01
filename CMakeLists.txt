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

find_package(Boost 1.39.0 QUIET REQUIRED COMPONENTS
  date_time
  filesystem
  log
  program_options
  system
  thread
  unit_test_framework
)
find_package(CUDA QUIET)
find_package(HDF5 QUIET REQUIRED)
find_package(Lua51 QUIET REQUIRED)
find_package(Luabind QUIET REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(NOT CMAKE_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "-fPIC -Wall")
  endif(NOT CMAKE_CXX_FLAGS)
endif()

if(CUDA_FOUND AND NOT HALMD_USE_STATIC_LIBS)
  set(WITH_CUDA TRUE)
  enable_language(CUDA)

  if(NOT CMAKE_CUDA_FLAGS)
    set(CMAKE_CUDA_FLAGS "-Xcompiler -fPIC -Xptxas -v --host-compilation=c -arch sm_12")
  endif(NOT CMAKE_CUDA_FLAGS)

  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}" CACHE STRING
    "Flags used by the compiler during all build types." FORCE)

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

set(HALMD_COMMON_LIBRARIES
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_LOG_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
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

else(NOT HALMD_DOC_ONLY)

project(HALMD NONE) # no compiler languages

add_subdirectory(doc)

endif(NOT HALMD_DOC_ONLY)

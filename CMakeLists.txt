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

set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/platform.cmake")

project(HALMD CXX)

if(HALMD_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS TRUE)
  set(HDF5_USE_STATIC_LIBS TRUE)
endif()

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

# Define macro LUABIND_DYNAMIC_LINK if linking to shared Luabind library.
#
# http://www.rasterbar.com/products/luabind/docs.html#building-and-testing
#
get_filename_component(suffix "${LUABIND_LIBRARY}" EXT)
if(suffix STREQUAL "${CMAKE_SHARED_LIBRARY_SUFFIX}")
  add_definitions(-DLUABIND_DYNAMIC_LINK)
endif()

if(CUDA_FOUND)
  enable_language(CUDA)

  set(HALMD_WITH_GPU TRUE)

  # FIXME deprecated, use HALMD_WITH_GPU instead
  set(WITH_CUDA TRUE)
  add_definitions(-DWITH_CUDA)

  try_compile(success "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/cmake/cuda_arch.cu" OUTPUT_VARIABLE result)
  if(result MATCHES "__CUDA_ARCH__([0-9]+)__")
    set(HALMD_GPU_ARCH "${CMAKE_MATCH_1}")
  else()
    message(SEND_ERROR "${result}")
  endif()

endif()

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

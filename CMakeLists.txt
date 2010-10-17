cmake_minimum_required(VERSION 2.8)

project(HALMD)

if (HALMD_USE_STATIC_LIBS)
  set(HALMD_USE_STATIC_LIBS TRUE CACHE INTERNAL
    "Build statically linked backend executables")
  set(Boost_USE_STATIC_LIBS TRUE)
  set(HDF5_USE_STATIC_LIBS TRUE)
  set(GSL_USE_STATIC_LIBS TRUE)
endif(HALMD_USE_STATIC_LIBS)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Boost 1.42.0 QUIET REQUIRED COMPONENTS
  date_time
  filesystem
  log
  program_options
  regex
  system
  thread
  unit_test_framework
  )
find_package(CUDA QUIET)
find_package(Git QUIET REQUIRED)
find_package(GSL QUIET REQUIRED)
find_package(HDF5 QUIET REQUIRED)
find_package(Lua51 QUIET REQUIRED)
find_package(Luabind QUIET REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Choose the type of build, options are: None Debug DeviceEmu Release RelWithDebInfo MinSizeRel." FORCE)

if(NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "-fPIC -Wall")
endif(NOT CMAKE_CXX_FLAGS)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING
  "Flags used by the compiler during all build types." FORCE)

set(CMAKE_CXX_FLAGS_DEVICEEMU "-O2 -g -D__DEVICE_EMULATION__" CACHE STRING
  "Flags used by the compiler during device emulation builds." FORCE)

set(CMAKE_C_FLAGS_DEVICEEMU "-O2 -g -D__DEVICE_EMULATION__" CACHE STRING
  "Flags used by the compiler during device emulation builds." FORCE)

set(CMAKE_EXE_LINKER_FLAGS_DEVICEEMU "" CACHE STRING
  "Flags used by the linker during device emulation builds." FORCE)

# Link shared modules with keyword 'nodelete', so the module does not get
# unloaded on dlclose. This works around a segmentation fault in the CUDA
# runtime library with device emulation.
set(CMAKE_MODULE_LINKER_FLAGS_DEVICEEMU "-z nodelete" CACHE STRING
  "Flags used by the linker during device emulation builds." FORCE)

set(CMAKE_SHARED_LINKER_FLAGS_DEVICEEMU "" CACHE STRING
  "Flags used by the shared libraries linker during device emulation builds." FORCE)

mark_as_advanced(
  CMAKE_CXX_FLAGS_DEVICEEMU
  CMAKE_C_FLAGS_DEVICEEMU
  CMAKE_EXE_LINKER_FLAGS_DEVICEEMU
  CMAKE_MODULE_LINKER_FLAGS_DEVICEEMU
  CMAKE_SHARED_LINKER_FLAGS_DEVICEEMU
  )

if(CMAKE_BUILD_TYPE STREQUAL "DeviceEmu")
  set(CUDA_EMULATION TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "DeviceEmu")

if(CUDA_EMULATION AND CUDA_RUNTIME_LIBRARY AND CUDA_INCLUDE_DIR AND NOT CUDA_FOUND)
  set(CUDA_LIBRARIES "${CUDA_RUNTIME_LIBRARY}")
  set(CUDA_FOUND "YES")
endif(CUDA_EMULATION AND CUDA_RUNTIME_LIBRARY AND CUDA_INCLUDE_DIR AND NOT CUDA_FOUND)

if(CUDA_FOUND)
  set(CMAKE_CUDA_FLAGS_DEVICEEMU "-O2 -g -deviceemu" CACHE STRING
   "Flags used by the compiler during device emulation builds." FORCE)

  mark_as_advanced(
    CMAKE_CUDA_FLAGS_DEVICEEMU
  )
endif(CUDA_FOUND)

if(CUDA_FOUND AND NOT HALMD_USE_STATIC_LIBS)
  set(WITH_CUDA TRUE)
  enable_language(CUDA)

  if(NOT CMAKE_CUDA_FLAGS)
    set(CMAKE_CUDA_FLAGS "-Xcompiler -fPIC -Xptxas -v --host-compilation=c")
  endif(NOT CMAKE_CUDA_FLAGS)

  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}" CACHE STRING
    "Flags used by the compiler during all build types." FORCE)

  add_definitions(-DWITH_CUDA)
else(CUDA_FOUND AND NOT HALMD_USE_STATIC_LIBS)
  set(WITH_CUDA FALSE)
endif(CUDA_FOUND AND NOT HALMD_USE_STATIC_LIBS)

git_repository_info(${HALMD_SOURCE_DIR} HALMD)

string(TOLOWER "${PROJECT_NAME}" PROGRAM_NAME)
set(PROGRAM_DESC "HAL’s MD package")
if(HALMD_GIT_VERSION)
  string(REGEX MATCH "v[0-9].*" PROGRAM_VERSION "${HALMD_GIT_VERSION}")
  if(NOT PROGRAM_VERSION)
    set(PROGRAM_VERSION "commit ${HALMD_GIT_VERSION}")
  endif(NOT PROGRAM_VERSION)
else(HALMD_GIT_VERSION)
  set(PROGRAM_VERSION "(devel)")
endif(HALMD_GIT_VERSION)
if(HALMD_GIT_AUTHOR_DATE)
  set(PROGRAM_DATE "${HALMD_GIT_AUTHOR_DATE}")
else(HALMD_GIT_AUTHOR_DATE)
  set(PROGRAM_DATE "unknown date")
endif(HALMD_GIT_AUTHOR_DATE)

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

set(HALMD_MODULES
  halmd_io
  halmd_io_profile
  halmd_io_profile_writers
  halmd_io_statevars
  halmd_io_statevars_writers
  halmd_io_trajectory
  halmd_io_trajectory_readers
  halmd_io_trajectory_writers
  halmd_mdsim
  halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_position
  halmd_mdsim_host_sample
  halmd_mdsim_host_sort
  halmd_mdsim_host_velocities
  halmd_mdsim_samples_host
  halmd_observables
  halmd_observables_host
  halmd_random
  halmd_random_host
  halmd_utility
  )
if(WITH_CUDA)
  set(HALMD_MODULES
    ${HALMD_MODULES}
    halmd_mdsim_gpu
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_position
    halmd_mdsim_gpu_sample
    halmd_mdsim_gpu_sort
    halmd_mdsim_gpu_velocities
    halmd_mdsim_samples_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
    )
endif(WITH_CUDA)

set(HALMD_COMMON_LIBRARIES
  -Wl,-whole-archive
  ${HALMD_MODULES}
  -Wl,-no-whole-archive
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${LUABIND_LIBRARIES}
  ${LUA_LIBRARIES}
  )

if(HALMD_USE_STATIC_LIBS)
  set(HALMD_COMMON_LIBRARIES
    ${HALMD_COMMON_LIBRARIES}
    -static-libgcc
    -Wl,-Bstatic
    )
endif(HALMD_USE_STATIC_LIBS)

set(HALMD_COMMON_LIBRARIES
  ${HALMD_COMMON_LIBRARIES}
  rt
  dl
  pthread
  z
  )

enable_testing()

include_directories(${HALMD_SOURCE_DIR}/include)
include_directories(${HALMD_SOURCE_DIR})
include_directories(${HALMD_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${LUA_INCLUDE_DIR})
include_directories(${LUABIND_INCLUDE_DIR})
if(WITH_CUDA)
  include_directories(${CUDA_INCLUDE_DIR})
endif(WITH_CUDA)

add_subdirectory(doc EXCLUDE_FROM_ALL)
add_subdirectory(halmd)
add_subdirectory(test)

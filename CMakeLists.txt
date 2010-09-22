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

set(PROGRAM_NAME "${PROJECT_NAME}")
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

enable_testing()

include_directories(${HALMD_SOURCE_DIR}/include)
include_directories(${HALMD_SOURCE_DIR}/src)
include_directories(${HALMD_BINARY_DIR}/src)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${LUA_INCLUDE_DIR})
include_directories(${LUABIND_INCLUDE_DIR})
if(WITH_CUDA)
  include_directories(${CUDA_INCLUDE_DIR})
endif(WITH_CUDA)

add_subdirectory(doc)
add_subdirectory(src)

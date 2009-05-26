cmake_minimum_required(VERSION 2.6)

project(ljgpu)

if (ljgpu_USE_STATIC_LIBS)
  set(ljgpu_USE_STATIC_LIBS TRUE CACHE INTERNAL
    "Build statically linked backend executables")
  set(Boost_USE_STATIC_LIBS TRUE)
  set(HDF5_USE_STATIC_LIBS TRUE)
  set(GSL_USE_STATIC_LIBS TRUE)
endif(ljgpu_USE_STATIC_LIBS)

set(CMAKE_MODULE_PATH ${ljgpu_SOURCE_DIR}/CMakeModules)

find_package(Boost QUIET REQUIRED COMPONENTS filesystem system date_time program_options unit_test_framework)
find_package(CUDA QUIET)
find_package(Git QUIET REQUIRED)
find_package(GSL QUIET REQUIRED)
find_package(HDF5 QUIET REQUIRED)

if(CUDA_FOUND AND NOT ljgpu_USE_STATIC_LIBS)
  set(WITH_CUDA TRUE)
  enable_language(CUDA)
  add_definitions(-DWITH_CUDA)
else(CUDA_FOUND AND NOT ljgpu_USE_STATIC_LIBS)
  set(WITH_CUDA FALSE)
endif(CUDA_FOUND AND NOT ljgpu_USE_STATIC_LIBS)

git_repository_info(${ljgpu_SOURCE_DIR} ljgpu)

set(PROGRAM_NAME "ljgpu")
if(ljgpu_GIT_VERSION)
  string(REGEX MATCH "v[0-9].*" PROGRAM_VERSION "${ljgpu_GIT_VERSION}")
  if(NOT PROGRAM_VERSION)
    set(PROGRAM_VERSION "(commit ${ljgpu_GIT_COMMIT})")
  endif(NOT PROGRAM_VERSION)
else(ljgpu_GIT_VERSION)
  set(PROGRAM_VERSION "(devel)")
endif(ljgpu_GIT_VERSION)
if(ljgpu_GIT_AUTHOR_DATE)
  set(PROGRAM_DATE "${ljgpu_GIT_AUTHOR_DATE}")
else(ljgpu_GIT_AUTHOR_DATE)
  set(PROGRAM_DATE "unknown date")
endif(ljgpu_GIT_AUTHOR_DATE)

enable_testing()

include_directories(${ljgpu_SOURCE_DIR}/include)
include_directories(${ljgpu_SOURCE_DIR}/src)
include_directories(${ljgpu_BINARY_DIR}/src)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
if(WITH_CUDA)
  include_directories(${CUDA_INCLUDE_DIR})
endif(WITH_CUDA)

add_subdirectory(src)

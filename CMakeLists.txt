cmake_minimum_required(VERSION 2.6)

project(ljgpu)

set(CMAKE_MODULE_PATH ${ljgpu_SOURCE_DIR}/CMakeModules)

find_package(Boost QUIET REQUIRED COMPONENTS filesystem date_time program_options)
find_package(CUDA QUIET)
find_package(Git QUIET REQUIRED)
find_package(GSL QUIET REQUIRED)
find_package(HDF5 QUIET REQUIRED)

if(CUDA_FOUND)
  set(WITH_CUDA TRUE)
  enable_language(CUDA)
  add_definitions(-DWITH_CUDA)
else(CUDA_FOUND)
  set(WITH_CUDA FALSE)
endif(CUDA_FOUND)

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

include_directories(${ljgpu_SOURCE_DIR}/include)
include_directories(${ljgpu_SOURCE_DIR}/src)
include_directories(${ljgpu_BINARY_DIR}/src)

add_subdirectory(src)

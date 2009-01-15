cmake_minimum_required(VERSION 2.6)

project(ljgpu)

set(CMAKE_MODULE_PATH ${ljgpu_SOURCE_DIR}/CMakeModules)

set(GIT_FIND_REQUIRED TRUE)
find_package(Git)

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


set(GSL_FIND_REQUIRED TRUE)
set(GSL_FIND_QUIETLY TRUE)
find_package(GSL)

set(Boost_FIND_REQUIRED TRUE)
set(Boost_FIND_QUIETLY TRUE)
find_package(Boost COMPONENTS filesystem date_time program_options)

enable_language(CUDA)

set(CUDA_FIND_REQUIRED TRUE)
set(CUDA_FIND_QUIETLY TRUE)
find_package(CUDA)

set(HDF5_FIND_REQUIRED TRUE)
set(HDF5_FIND_QUIETLY TRUE)
find_package(HDF5)

include_directories(${ljgpu_SOURCE_DIR}/include)
include_directories(${ljgpu_SOURCE_DIR}/src)
include_directories(${ljgpu_BINARY_DIR}/src)

add_subdirectory(src)
add_subdirectory(test)

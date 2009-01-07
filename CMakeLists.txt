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

# Add program variant compile flags to cmake GUI
set(ljgpu_VARIANT_CUDA TRUE CACHE BOOL
  "Use CUDA")
set(ljgpu_VARIANT_CELL FALSE CACHE BOOL
  "Use fixed-size cell lists")
set(ljgpu_VARIANT_CELL_SUMMATION_ORDER FALSE CACHE BOOL
  "Use opposite cell summation order")
set(ljgpu_VARIANT_DSFUN FALSE CACHE BOOL
  "Use double-single precision functions")
set(ljgpu_VARIANT_NEIGHBOUR TRUE CACHE BOOL
  "Use Verlet neighbour lists")
set(ljgpu_VARIANT_HILBERT TRUE CACHE BOOL
  "Use Hilbert space-filling curve particle ordering")
set(ljgpu_VARIANT_C2POT FALSE CACHE BOOL
  "Use C2 potential smoothing")

if(ljgpu_VARIANT_CUDA)
  add_definitions(-DUSE_CUDA)
  if(ljgpu_VARIANT_NEIGHBOUR)
    add_definitions(-DUSE_NEIGHBOUR)
    if(ljgpu_VARIANT_HILBERT)
      add_definitions(-DUSE_HILBERT_ORDER)
    endif(ljgpu_VARIANT_HILBERT)
  else(ljgpu_VARIANT_NEIGHBOUR)
    if(ljgpu_VARIANT_CELL)
      add_definitions(-DUSE_CELL)
      if(ljgpu_VARIANT_CELL_SUMMATION_ORDER)
	add_definitions(-DUSE_CELL_SUMMATION_ORDER)
      endif(ljgpu_VARIANT_CELL_SUMMATION_ORDER)
      if(ljgpu_VARIANT_DSFUN)
	add_definitions(-DUSE_DSFUN)
      endif(ljgpu_VARIANT_DSFUN)
    endif(ljgpu_VARIANT_CELL)
  endif(ljgpu_VARIANT_NEIGHBOUR)
endif(ljgpu_VARIANT_CUDA)
if(ljgpu_VARIANT_C2POT)
  add_definitions(-DUSE_POTENTIAL_SMOOTHING)
endif(ljgpu_VARIANT_C2POT)

include_directories(${ljgpu_SOURCE_DIR}/include)
include_directories(${ljgpu_SOURCE_DIR}/src)
include_directories(${ljgpu_BINARY_DIR}/src)

add_subdirectory(src)
add_subdirectory(test)

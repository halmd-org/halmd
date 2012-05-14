cmake_minimum_required(VERSION 2.8.8)

if(NOT HALMD_DOC_ONLY)
  set(HALMD_DOC_ONLY FALSE)
else()
  set(HALMD_DOC_ONLY TRUE)
endif()
set(HALMD_DOC_ONLY ${HALMD_DOC_ONLY} CACHE INTERNAL "Build HALMD documentation only" FORCE)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# If the subdirectory cmake/ contains a version.cmake.in file, we are in a git
# repository, then extract the git commit from the repository to generate
# CMake files with version information and rules for crafting a HALMD archive.
#
# If the subdirectory cmake/ does not contain version.cmake.in, we are in a
# HALMD archive, then include pre-build CMake file with version information.
#
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/version.cmake.in")
  find_package(Git QUIET REQUIRED)
  find_package(GitRepository QUIET REQUIRED)

  git_repository("${CMAKE_SOURCE_DIR}" HALMD)

  configure_file(cmake/version.cmake.in cmake/version.cmake @ONLY)
  include("${CMAKE_BINARY_DIR}/cmake/version.cmake")
  include("${CMAKE_SOURCE_DIR}/cmake/archive.cmake")
else()
  include("${CMAKE_SOURCE_DIR}/cmake/version.cmake")
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

# Require thread-safe Boost libraries for thread-safe Boost Log
set(Boost_USE_MULTITHREADED TRUE)

find_package(Boost 1.43.0 QUIET REQUIRED COMPONENTS
  date_time
  filesystem
  program_options
  random
  system
  thread
  unit_test_framework
)
find_package(CUDA QUIET)
find_package(HDF5 QUIET REQUIRED)
find_package(Lua QUIET REQUIRED)

if(CUDA_FOUND)
  enable_language(CUDA)

  set(HALMD_WITH_GPU TRUE)

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

if(HALMD_WITH_GPU)
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

endif(HALMD_WITH_GPU)

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
  boost_log
  luabind
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_RANDOM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  ${HDF5_CPP_LIBRARY}
  ${HDF5_LIBRARY}
  ${LUA_LIBRARIES}
)
if(HALMD_WITH_GPU)
  list(APPEND HALMD_COMMON_LIBRARIES
    ${CUDA_LIBRARIES}
  )
endif(HALMD_WITH_GPU)
list(APPEND HALMD_COMMON_LIBRARIES
  rt
  dl
  pthread
  z
)

enable_testing()
include(CTest)

include_directories(${HALMD_SOURCE_DIR}/libs/boost-log)
include_directories(${HALMD_SOURCE_DIR}/libs/cuda-wrapper)
include_directories(${HALMD_SOURCE_DIR}/libs/h5xx)
include_directories(${HALMD_SOURCE_DIR}/libs/luabind)
include_directories(${HALMD_SOURCE_DIR})
include_directories(${HALMD_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${LUA_INCLUDE_DIR})
if(HALMD_WITH_GPU)
  include_directories(${CUDA_INCLUDE_DIR})
endif(HALMD_WITH_GPU)

# If the subdirectory doc/ contains a CMakeLists.txt, we are in a git repository,
# then include the CMake rules in doc/ to extract and generate documentation.
# Note that the targets halmd_doc_html, halmd_doc_pdf and halmd_doc_man are
# not included in the default targets, so to generate documentation, one
# has to explicitly invoke the desired target with make halmd_doc_<format>.
#
# If the subdirectory doc/ does not contain a CMakeLists.txt, we are in a
# HALMD archive, then add a rule to install the pre-build documentation.
#
if(EXISTS "${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt")
  add_subdirectory(doc)
else()
  install(DIRECTORY "${CMAKE_SOURCE_DIR}/doc/"
    DESTINATION "share/doc/halmd"
    PATTERN "man" EXCLUDE
  )
  install(DIRECTORY "${CMAKE_SOURCE_DIR}/doc/man/"
    DESTINATION "share/man/man1"
  )
endif()

# Install simulation and plotting examples to <prefix>/share/doc/halmd.
# This is common practice in GNU/Linux distributions, and encourages
# users to make a copy of an example and modify it for their needs.
#
install(DIRECTORY ${CMAKE_SOURCE_DIR}/examples
  DESTINATION share/doc/halmd
  PATTERN packages.mk EXCLUDE
)

add_subdirectory(libs/boost-log/libs/log/src)
add_subdirectory(libs/luabind/src)
add_subdirectory(halmd)
add_subdirectory(lua)
add_subdirectory(test)
add_subdirectory(libs/h5xx/test)

else(NOT HALMD_DOC_ONLY)

project(HALMD NONE) # no compiler languages

add_subdirectory(doc)

endif(NOT HALMD_DOC_ONLY)

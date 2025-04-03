cmake_minimum_required(VERSION 3.8)
project(HALMD NONE)

set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/platform.cmake")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set cmake policies to avoid warnings in newer cmake versions.
if(POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif()
if(POLICY CMP0026)
  cmake_policy(SET CMP0026 NEW)
endif()
if(POLICY CMP0058)
  cmake_policy(SET CMP0058 NEW)
endif()

# If the subdirectory cmake/ contains a version.cmake.in file, we are in a git
# repository, then extract the git commit from the repository to generate
# CMake files with version information and rules for crafting a HALMD archive.
#
# If the subdirectory cmake/ does not contain version.cmake.in, we are in a
# HALMD archive, then include pre-build CMake file with version information.
#
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/version.cmake.in")
  set(HALMD_USE_GIT TRUE)
else()
  set(HALMD_USE_GIT FALSE)
endif()

if(HALMD_USE_GIT)
  find_package(Git QUIET REQUIRED)
  find_package(GitRepository QUIET REQUIRED)

  git_repository("${CMAKE_SOURCE_DIR}" HALMD)

  configure_file(cmake/version.cmake.in cmake/version.cmake @ONLY)
  include("${CMAKE_BINARY_DIR}/cmake/version.cmake")

  # Target to generate HALMD archive
  include("${CMAKE_SOURCE_DIR}/cmake/archive.cmake")

else()
  include("${CMAKE_SOURCE_DIR}/cmake/version.cmake")
endif()

# Output source version for Dashboard test reports
message(STATUS "Building HALMD ${PROGRAM_VERSION}")

# Set HALMD potential list needed for building everything or only the docs
include(${CMAKE_SOURCE_DIR}/cmake/potentials.cmake)

if(HALMD_USE_GIT)
  if(HALMD_DOC_ONLY)
    set(HALMD_DOC_ONLY TRUE)
  else()
    set(HALMD_DOC_ONLY FALSE)
  endif()
  set(HALMD_DOC_ONLY ${HALMD_DOC_ONLY} CACHE INTERNAL "Build HALMD documentation only" FORCE)
endif()

if(NOT HALMD_DOC_ONLY OR NOT HALMD_USE_GIT)
  # By default, enable GPU support only if CUDA is available.
  # If HALMD_WITH_GPU is explicitly set to TRUE, require CUDA.
  include(CheckLanguage)
  check_language(CUDA)

  if(CMAKE_CUDA_COMPILER)
    enable_language(CXX CUDA)

    set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "Use the static version of the CUDA runtime library if available")

    if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 9.0)
      message(SEND_ERROR "Minium required CUDA version is 9.0 (found ${CMAKE_CUDA_COMPILER_VERSION})")
    endif()

    # warn user about static linking when using CUDA 9.2 because of a bug
    if(${CUDA_USE_STATIC_CUDA_RUNTIME} AND ${CMAKE_CUDA_COMPILER_VERSION} VERSION_EQUAL 9.2)
      message(WARNING "Static linking of the CUDA runtime library may lead to kernel launch failures")
    endif()

    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD 14)

    set(CMAKE_CUDA_ARCHITECTURES OFF)

    find_library(CUDA_LIBRARY NAMES cuda HINTS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    if(CUDA_USE_STATIC_CUDA_RUNTIME)
      find_library(CUDART_LIBRARY cudart_static HINTS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    else()
      find_library(CUDART_LIBRARY cudart HINTS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    endif()

    if(NOT CUDA_LIBRARY)
      message(SEND_ERROR "CUDA driver library could not be found")
    elseif(NOT CUDART_LIBRARY)
      message(SEND_ERROR "CUDA runtime library could not be found")
    endif()

    # unset CUDA implicit libraries so they don't overwrite the C++ libraries
    unset(CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES)
    unset(CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES)
    # we manually choose between cudart static/dynamic
    set(CMAKE_CUDA_RUNTIME_LIBRARY None)

    try_compile(success "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/cmake/cuda_arch.cu" OUTPUT_VARIABLE result)
    if(result MATCHES "__CUDA_ARCH__([0-9]+)__")
      set(HALMD_GPU_ARCH "${CMAKE_MATCH_1}")
    else()
      message(SEND_ERROR "${result}")
    endif()

    set(HALMD_WITH_GPU TRUE)
  else()
    if(HALMD_WITH_GPU)
      message(SEND_ERROR "CUDA not found")
    endif()

    enable_language(CXX)

    set(HALMD_WITH_GPU FALSE)
  endif()

  set(HALMD_USE_STATIC_LIBS FALSE CACHE BOOL
      "Use static linkage for Boost, HDF5, and Lua libraries"
  )

  if(HALMD_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS TRUE)
    set(HDF5_USE_STATIC_LIBS TRUE)
    set(LUA_USE_STATIC_LIBS TRUE)

    message(STATUS "For static linking of the CUDA runtime library specifically set CUDA_USE_STATIC_CUDA_RUNTIME to ON")
  endif()

  # define BOOST_LOG_DYN_LINK in case of dynamic linkage
  if(NOT Boost_USE_STATIC_LIBS)
     add_definitions(-DBOOST_LOG_DYN_LINK)
  endif()

  # Require thread-safe Boost libraries
  set(Boost_USE_MULTITHREADED TRUE)

  find_package(Boost 1.66.0 QUIET REQUIRED COMPONENTS
    date_time
    filesystem
    log
    program_options
    random
    system
    thread
    unit_test_framework
  )
  set(Boost_VERSION ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION})

  find_package(HDF5 QUIET REQUIRED COMPONENTS CXX)
  find_package(LuaLibs QUIET REQUIRED)

  # detect HDF5 version manually for cmake < 3.3
  if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS 3.3)
    find_path(H5PUBCONF_PATH H5pubconf.h PATHS ${HDF5_INCLUDE_DIRS} NO_DEFAULT_PATH)
    file(STRINGS "${H5PUBCONF_PATH}/H5pubconf.h" HDF5_VERSION LIMIT_COUNT 1 REGEX "#define H5_VERSION \"(.*)\"")
    string(REGEX REPLACE "#define H5_VERSION \"(.*)\"" "\\1" HDF5_VERSION "${HDF5_VERSION}")
  endif()

  # be verbose about library versions (but suppress all the gory details given by find_package)
  message(STATUS "Boost version: ${Boost_VERSION}")
  message(STATUS "Lua library version: ${LUA_VERSION_STRING}")
  message(STATUS "HDF5 library version: ${HDF5_VERSION}")

  # Set HALMD build variant flags
  include(${CMAKE_SOURCE_DIR}/cmake/variant.cmake)

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
    luaponte
    ${Boost_DATE_TIME_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_RANDOM_LIBRARY}
    ${Boost_LOG_LIBRARY} # must precede Boost.Filesystem for static linkage
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_THREAD_LIBRARY}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${LUA_LIBRARIES}
  )
  if(HALMD_WITH_GPU)
    list(APPEND HALMD_COMMON_LIBRARIES
      ${CUDART_LIBRARY}
      ${CUDA_LIBRARY}
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

  include_directories(${HALMD_SOURCE_DIR}/libs/h5xx)
  include_directories(${HALMD_SOURCE_DIR}/libs/luaponte)
  include_directories(${HALMD_SOURCE_DIR})
  include_directories(${HALMD_BINARY_DIR})
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
  include_directories(SYSTEM ${LUA_INCLUDE_DIR})
  if(HALMD_WITH_GPU)
    include_directories(${HALMD_SOURCE_DIR}/libs/cub)
    include_directories(${HALMD_SOURCE_DIR}/libs/cuda-wrapper)
    include_directories(SYSTEM ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

    if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0)
      add_definitions(-DTHRUST_IGNORE_CUB_VERSION_CHECK)
    endif()
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
  if(HALMD_USE_GIT)
    add_subdirectory(doc)
  elseif(EXISTS "${CMAKE_SOURCE_DIR}/doc/")
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

  add_subdirectory(libs/luaponte/src)
  add_subdirectory(halmd)
  add_subdirectory(lua)
  add_subdirectory(test)
  if(HALMD_WITH_GPU)
    add_subdirectory(libs/cuda-wrapper/test)
  endif()
  add_subdirectory(libs/h5xx/test)
  add_subdirectory(libs/luaponte/test)

else()

  # Set HALMD build variant flags
  include(${CMAKE_SOURCE_DIR}/cmake/variant.cmake)

  add_subdirectory(doc)
  add_subdirectory(lua)

endif()

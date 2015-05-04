cmake_minimum_required(VERSION 2.8.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

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

  set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/platform.cmake")

  project(HALMD CXX)

  # By default, enable GPU support only if CUDA is available.
  # If HALMD_WITH_GPU is explicitly set to TRUE, require CUDA.
  if(NOT DEFINED HALMD_WITH_GPU)
    find_package(CUDA QUIET COMPONENTS cudart)
  elseif(HALMD_WITH_GPU)
    find_package(CUDA QUIET REQUIRED COMPONENTS cudart)
  endif()

  if(CUDA_FOUND)
    enable_language(CUDA)
    if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 5.0)
      message(SEND_ERROR "Minium required CUDA version is 5.0 (found ${CMAKE_CUDA_COMPILER_VERSION})")
    endif()

    set(HALMD_WITH_GPU TRUE)

    try_compile(success "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/cmake/cuda_arch.cu" OUTPUT_VARIABLE result)
    if(result MATCHES "__CUDA_ARCH__([0-9]+)__")
      set(HALMD_GPU_ARCH "${CMAKE_MATCH_1}")
    else()
      message(SEND_ERROR "${result}")
    endif()
  else()
    set(HALMD_WITH_GPU FALSE)
  endif()

  set(HALMD_USE_STATIC_LIBS FALSE CACHE BOOL
      "Use static linkage for Boost, HDF5, and Lua libraries"
  )

  if(HALMD_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS TRUE)
    set(HDF5_USE_STATIC_LIBS TRUE)
    set(LUA_USE_STATIC_LIBS TRUE)
  endif()

  # define BOOST_LOG_DYN_LINK in case of dynamic linkage
  if(NOT Boost_USE_STATIC_LIBS)
     add_definitions(-DBOOST_LOG_DYN_LINK)
  endif()

  # Require thread-safe Boost libraries
  set(Boost_USE_MULTITHREADED TRUE)

  find_package(Boost 1.55.0 QUIET REQUIRED COMPONENTS
    date_time
    filesystem
    log
    program_options
    random
    system
    thread
    unit_test_framework
  )
  find_package(HDF5 QUIET REQUIRED)
  find_package(LuaLibs QUIET REQUIRED)

  # be verbose about library versions (but suppress all the gory details given by find_package)
  message(STATUS "Boost version: ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
  message(STATUS "Lua library version: ${LUA_VERSION_STRING}")
  message(STATUS "HDF5 library version: ${HDF5_VERSION_STRING}")

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

  include_directories(${HALMD_SOURCE_DIR}/libs/cuda-wrapper)
  include_directories(${HALMD_SOURCE_DIR}/libs/h5xx)
  include_directories(${HALMD_SOURCE_DIR}/libs/luaponte)
  include_directories(${HALMD_SOURCE_DIR})
  include_directories(${HALMD_BINARY_DIR})
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  include_directories(SYSTEM ${HDF5_INCLUDE_DIR})
  include_directories(SYSTEM ${LUA_INCLUDE_DIR})
  if(HALMD_WITH_GPU)
    include_directories(SYSTEM ${CUDA_INCLUDE_DIR})
    include_directories(${HALMD_SOURCE_DIR}/libs/cub)
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

  project(HALMD NONE) # no compiler languages

  add_subdirectory(doc)
  add_subdirectory(lua)

endif()

# - Find the Luabind headers and library
#
# This module defines
#  LUABIND_FOUND
#  LUABIND_INCLUDE_DIR
#  LUABIND_LIBRARY

#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
# Copyright 2010-2011 Peter Colberg
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file COPYING-CMAKE-SCRIPTS for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_path(LUABIND_INCLUDE_DIR luabind/luabind.hpp
  HINTS
  $ENV{LUABIND_DIR}
  PATH_SUFFIXES include
)

if(LUABIND_USE_STATIC_LIBS)
  set(_LUABIND_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
endif()

if(LUABIND_USE_DEBUG_LIBS)
  set(_LUABIND_LIBRARY_NAME luabindd)
else()
  set(_LUABIND_LIBRARY_NAME luabind)
endif()

find_library(LUABIND_LIBRARY NAMES ${_LUABIND_LIBRARY_NAME}
  HINTS
  $ENV{LUABIND_DIR}
  PATH_SUFFIXES lib64 lib
)

if(LUABIND_USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_LUABIND_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Luabind DEFAULT_MSG
  LUABIND_INCLUDE_DIR
  LUABIND_LIBRARY
)

mark_as_advanced(
  LUABIND_INCLUDE_DIR
  LUABIND_LIBRARY
)

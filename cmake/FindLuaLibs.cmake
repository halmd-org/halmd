# Locate Lua library
# This module defines
#  LUALIBS_FOUND
#  LUA_LIBRARIES
#  LUA_INCLUDE_DIR
#  LUA_VERSION_STRING
#
# Note that the expected include convention is
#  #include "lua.h"
# and not
#  #include <lua/lua.h>
# This is because, the lua location is not standardized and may exist
# in locations other than lua/

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
# Copyright 2011 Peter Colberg
# Copyright 2013 Felix Höfling
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_path(LUA_INCLUDE_DIR lua.h
  HINTS
    ENV LUA_DIR
  PATH_SUFFIXES include/luajit-2.0 include/lua52 include/lua5.2 include/lua-5.2 include/lua51 include/lua5.1 include/lua-5.1 include/lua include
  PATHS
  ~/Library/Frameworks
  /Library/Frameworks
  /sw # Fink
  /opt/local # DarwinPorts
  /opt/csw # Blastwave
  /opt
)

if(LUA_USE_STATIC_LIBS)
  set( _LUA_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
endif()

find_library(LUA_LIBRARY
  NAMES luajit-5.2 luajit-5.1 lua lua52 lua5.2 lua-5.2 lua51 lua5.1 lua-5.1
  HINTS
    ENV LUA_DIR
  PATH_SUFFIXES lib64 lib
  PATHS
  ~/Library/Frameworks
  /Library/Frameworks
  /sw
  /opt/local
  /opt/csw
  /opt
)

if(LUA_USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_LUA_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

if(LUA_LIBRARY)
  # include the math library for Unix
  if(UNIX AND NOT APPLE AND NOT BEOS)
    find_library(LUA_MATH_LIBRARY m)
    set( LUA_LIBRARIES "${LUA_LIBRARY};${LUA_MATH_LIBRARY}" CACHE STRING "Lua Libraries")
  # For Windows and Mac, don't need to explicitly include the math library
  else()
    set( LUA_LIBRARIES "${LUA_LIBRARY}" CACHE STRING "Lua Libraries")
  endif()
endif()

if(LUA_INCLUDE_DIR AND EXISTS "${LUA_INCLUDE_DIR}/lua.h")
  file(STRINGS "${LUA_INCLUDE_DIR}/lua.h" lua_version_str REGEX "^#define[ \t]+LUA_RELEASE[ \t]+\"Lua .+\"")

  string(REGEX REPLACE "^#define[ \t]+LUA_RELEASE[ \t]+\"Lua ([^\"]+)\".*" "\\1" LUA_VERSION_STRING "${lua_version_str}")
  unset(lua_version_str)

  # the above does not work for Lua 5.2
  if(LUA_VERSION_STRING STREQUAL "")
    file(STRINGS "${LUA_INCLUDE_DIR}/lua.h" lua_version_str REGEX "^#define[ \t]+LUA_VERSION_MAJOR[ \t]+\".+\"")
    string(REGEX REPLACE "^#define[ \t]+LUA_VERSION_MAJOR[ \t]+\"([^\"]+)\".*" "\\1" lua_version_major "${lua_version_str}")

    file(STRINGS "${LUA_INCLUDE_DIR}/lua.h" lua_version_str REGEX "^#define[ \t]+LUA_VERSION_MINOR[ \t]+\".+\"")
    string(REGEX REPLACE "^#define[ \t]+LUA_VERSION_MINOR[ \t]+\"([^\"]+)\".*" "\\1" lua_version_minor "${lua_version_str}")

    file(STRINGS "${LUA_INCLUDE_DIR}/lua.h" lua_version_str REGEX "^#define[ \t]+LUA_VERSION_RELEASE[ \t]+\".+\"")
    string(REGEX REPLACE "^#define[ \t]+LUA_VERSION_RELEASE[ \t]+\"([^\"]+)\".*" "\\1" lua_version_release "${lua_version_str}")

    set(LUA_VERSION_STRING "${lua_version_major}.${lua_version_minor}.${lua_version_release}")
    unset(lua_version_major)
    unset(lua_version_minor)
    unset(lua_version_release)
  endif()
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LUALIBS_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args(LuaLibs
                                  REQUIRED_VARS LUA_LIBRARIES LUA_INCLUDE_DIR
                                  VERSION_VAR LUA_VERSION_STRING)

mark_as_advanced(LUA_INCLUDE_DIR LUA_LIBRARIES LUA_LIBRARY LUA_MATH_LIBRARY)

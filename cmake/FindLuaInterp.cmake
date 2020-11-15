# Locate Lua interpreter
# This module defines
#  LUAINTERP_FOUND
#  LUA_EXECUTABLE
#

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
# Copyright 2011 Peter Colberg
# Copyright 2018 Felix Höfling
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

if(${CMAKE_VERSION} VERSION_LESS "3.4.0")
    message("For a consistent Lua configuration, please consider switching to CMake ≥ 3.4.0")
    find_program(LUA_EXECUTABLE
      NAMES luajit lua lua52 lua5.2 lua-5.2 lua51 lua5.1 lua-5.1
      # NAMES_PER_DIR   # not known to CMake < 3.4
      HINTS
        ENV LUA_DIR
      PATH_SUFFIXES bin
      PATHS
      ~/Library/Frameworks
      /Library/Frameworks
      /sw
      /opt/local
      /opt/csw
      /opt
    )
else()
    find_program(LUA_EXECUTABLE
      NAMES luajit lua lua52 lua5.2 lua-5.2 lua51 lua5.1 lua-5.1
      NAMES_PER_DIR
      HINTS
        ENV LUA_DIR
      PATH_SUFFIXES bin
      PATHS
      ~/Library/Frameworks
      /Library/Frameworks
      /sw
      /opt/local
      /opt/csw
      /opt
    )
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LUAINTERP_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args(LuaInterp DEFAULT_MSG LUA_EXECUTABLE)

mark_as_advanced(
  LUA_EXECUTABLE
)

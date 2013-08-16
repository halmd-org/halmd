# Locate Lua interpreter
# This module defines
#  LUAINTERP_FOUND
#  LUA_EXECUTABLE
#

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
# Copyright 2011 Peter Colberg
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

find_program(LUA_EXECUTABLE
  NAMES luajit lua lua52 lua5.2 lua-5.2 lua51 lua5.1 lua-5.1
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

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LUAINTERP_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args(LuaInterp DEFAULT_MSG LUA_EXECUTABLE)

mark_as_advanced(
  LUA_EXECUTABLE
)

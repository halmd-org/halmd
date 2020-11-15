# - This module looks for Sphinx
# Find the Sphinx documentation generator
#
# This modules defines
#  SPHINX_EXECUTABLE
#  SPHINX_FOUND

#=============================================================================
# Copyright 2016      Daniel Kirchner
# Copyright 2009-2011 Peter Colberg
# Copyright 2002-2009 Kitware, Inc.
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

find_program(SPHINX_EXECUTABLE NAMES sphinx-build
  HINTS
  $ENV{SPHINX_DIR}
  PATH_SUFFIXES bin
  DOC "Sphinx documentation generator"
)

if(SPHINX_EXECUTABLE)
   execute_process (COMMAND "${SPHINX_EXECUTABLE}" --version
                    OUTPUT_VARIABLE _SPHINX_VERSION
                    ERROR_VARIABLE  _SPHINX_VERSION
   )

   if (_SPHINX_VERSION MATCHES "sphinx-build\\)? ([0-9]+(\\.[0-9]+)+)")
     set (SPHINX_VERSION_STRING "${CMAKE_MATCH_1}")
  endif()
endif ()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx
  REQUIRED_VARS SPHINX_EXECUTABLE
  VERSION_VAR SPHINX_VERSION_STRING
)

mark_as_advanced(
  SPHINX_EXECUTABLE
)

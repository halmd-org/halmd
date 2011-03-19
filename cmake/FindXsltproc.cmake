# - This module looks for xsltproc
# Find the command line XSLT processor
#
# This modules defines
#  XSLTPROC_FOUND
#  XSLTPROC_EXECUTABLE

#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
# Copyright 2011 Peter Colberg
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

find_program(XSLTPROC_EXECUTABLE NAMES xsltproc
  HINTS
  $ENV{XSLTPROC_DIR}
  PATH_SUFFIXES bin
  DOC "XSLT processor"
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Xsltproc DEFAULT_MSG
  XSLTPROC_EXECUTABLE
)

mark_as_advanced(
  XSLTPROC_EXECUTABLE
)

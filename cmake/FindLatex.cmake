# - Find Latex
#
# This module defines
#  LATEX_FOUND
#  LATEX_COMPILER
#  PDFLATEX_COMPILER
#  DVIPNG_CONVERTER

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

find_program(LATEX_COMPILER NAMES latex
  HINTS
  $ENV{LATEX_DIR}
  PATH_SUFFIXES bin
  DOC "LaTeX compiler"
)

find_program(PDFLATEX_COMPILER NAMES pdflatex
  HINTS
  $ENV{LATEX_DIR}
  PATH_SUFFIXES bin
  DOC "pdfLaTeX compiler"
)

find_program(DVIPNG_CONVERTER NAMES dvipng
  HINTS
  $ENV{LATEX_DIR}
  PATH_SUFFIXES bin
  DOC "DVI-to-PNG converter"
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Latex DEFAULT_MSG
  LATEX_COMPILER
  PDFLATEX_COMPILER
  DVIPNG_CONVERTER
)

MARK_AS_ADVANCED(
  LATEX_COMPILER
  PDFLATEX_COMPILER
  DVIPNG_CONVERTER
)

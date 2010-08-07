# - This module looks for Sphinx
# Sphinx is a documentation generator.  Please see
# http://sphinx.pocoo.org/
#
# This modules defines the following variables:
#
#   SPHINX_EXECUTABLE     = The path to the sphinx-build command.
#   SPHINX_FOUND          = Was Sphinx found or not?
#

FIND_PROGRAM(SPHINX_EXECUTABLE
  NAMES sphinx-build
  DOC "Sphinx documentation generator (http://sphinx.pocoo.org/)"
)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)

MARK_AS_ADVANCED(
  SPHINX_EXECUTABLE
  )

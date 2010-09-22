# - Find the Luabind headers and library
#
# This module defines
#  LUABIND_INCLUDE_DIR
#  LUABIND_LIBRARIES
#  LUABIND_FOUND

find_path(LUABIND_INCLUDE_DIR
  luabind/luabind.hpp
)

find_library(LUABIND_LIBRARY
  NAMES luabind
)

find_library(LUABIND_LIBRARY_DEBUG
  NAMES luabindd
)

if(LUABIND_LIBRARY AND LUABIND_INCLUDE_DIR)
  set(LUABIND_LIBRARIES ${LUABIND_LIBRARY})
endif(LUABIND_LIBRARY AND LUABIND_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)

# Handle the REQUIRED and QUIET argument to FIND_PACKAGE()
# and set the LUABIND_FOUND variable.
find_package_handle_standard_args(Luabind DEFAULT_MSG
  LUABIND_LIBRARIES
  LUABIND_INCLUDE_DIR
)

mark_as_advanced(
  LUABIND_INCLUDE_DIR
  LUABIND_LIBRARY
  LUABIND_LIBRARY_DEBUG
)

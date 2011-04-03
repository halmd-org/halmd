##
# Add HALMD module target
#
# This macro substitutes the add_library command, and in addition
# adds the target name to the global property HALMD_MODULES.
#
# Base classes must be ordered *before* derived classes, otherwise
# Luabind will throw an assertion error or cause a segmentation fault.
#
macro(halmd_add_module MODULE)
  set_property(GLOBAL APPEND PROPERTY HALMD_MODULES ${MODULE})
  add_library(${MODULE} ${ARGN})
endmacro()

if(HALMD_WITH_GPU)
  add_definitions(-DHALMD_WITH_GPU)
endif(HALMD_WITH_GPU)

#
# The following option only works on x86-64 by default, which always uses the
# SSE instruction set for floating-point math. On i386, the x87 floating-point
# unit provides 80-bit extended double precision math internally, causing
# excess precision even if values are stored in single-precision.
#
# See the -mfpmath option in the gcc manpage for details, and
#
# Deterministic cross-platform floating point arithmetics
# http://www.christian-seiler.de/projekte/fpmath/
#
set(HALMD_VARIANT_HOST_SINGLE_PRECISION FALSE CACHE BOOL
  "Use single-precision math in host implementation (requires SSE)")
if(HALMD_VARIANT_HOST_SINGLE_PRECISION)
  add_definitions(-DUSE_HOST_SINGLE_PRECISION)
endif(HALMD_VARIANT_HOST_SINGLE_PRECISION)

if(HALMD_WITH_GPU OR HALMD_DOC_ONLY)
  set(HALMD_VARIANT_GPU_SINGLE_PRECISION FALSE CACHE BOOL
          "Enable single-precision math in gpu implementation")
  if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
    add_definitions(-DUSE_GPU_SINGLE_PRECISION)
  endif(HALMD_VARIANT_GPU_SINGLE_PRECISION)

  set(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION TRUE CACHE BOOL
          "Enable double-single-precision math in gpu implementation")
  if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    add_definitions(-DUSE_GPU_DOUBLE_SINGLE_PRECISION)
  endif(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)

  if(NOT HALMD_VARIANT_GPU_SINGLE_PRECISION AND NOT HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    message(SEND_ERROR "Either HALMD_VARIANT_GPU_SINGLE_PRECISION or HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION has to be set.")
  endif()
endif(HALMD_WITH_GPU OR HALMD_DOC_ONLY)

# cmake variables for string substitution in Lua files and documentation
if(HALMD_VARIANT_HOST_SINGLE_PRECISION)
  set(HALMD_HOST_FLOAT_TYPE "float")
  set(HALMD_HOST_PRECISION "single")
else()
  set(HALMD_HOST_FLOAT_TYPE "double")
  set(HALMD_HOST_PRECISION "double")
endif()

if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
  set(HALMD_DEFAULT_GPU_PRECISION "double-single")
else()
  set(HALMD_DEFAULT_GPU_PRECISION "single")
endif()

message(STATUS "Floating-point precision of host backend: ${HALMD_HOST_PRECISION}")
message(STATUS "Default floating-point precision of GPU backend: ${HALMD_DEFAULT_GPU_PRECISION}")

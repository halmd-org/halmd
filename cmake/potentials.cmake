if(NOT DEFINED HALMD_POTENTIALS)
  # find all potential modules
  file(GLOB potential_files "lua/halmd/mdsim/potentials/*.lua.in")
  foreach(file ${potential_files})
    # extract potential name from file name
    string(REGEX REPLACE ".*/(.*)\\.lua\\.in$" "\\1" result ${file})
    if(NOT ${result} MATCHES "init")
      list(APPEND HALMD_POTENTIALS ${result})
    endif()
    unset(result)
  endforeach(file)
endif()

# define cached variable, enable all potentials by default
set(HALMD_POTENTIALS ${HALMD_POTENTIALS} CACHE STRING
  "List of enabled potentials")

# define variables HALMD_WITH_*
foreach(potential ${HALMD_POTENTIALS})
  list(FIND "${HALMD_POTENTIALS}" ${potential} HALMD_WITH_${potential})
endforeach()

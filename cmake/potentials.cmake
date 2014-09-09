if(NOT DEFINED HALMD_EXTERNAL_POTENTIALS)
  # find all external potential modules
  file(GLOB potential_files "lua/halmd/mdsim/potentials/external/*.lua.in")
  foreach(file ${potential_files})
    # extract potential name from file name
    string(REGEX REPLACE ".*/external/(.*)\\.lua\\.in$" "\\1" result ${file})
    if(NOT ${result} MATCHES "init")
      list(APPEND HALMD_EXTERNAL_POTENTIALS ${result})
    endif()
    unset(result)
  endforeach(file)
endif()

if(NOT DEFINED HALMD_PAIR_POTENTIALS)
  # find all pair potential modules
  file(GLOB potential_files "lua/halmd/mdsim/potentials/pair/*.lua.in")
  foreach(file ${potential_files})
    # extract potential name from file name
    string(REGEX REPLACE ".*/pair/(.*)\\.lua\\.in$" "\\1" result ${file})
    if(NOT ${result} MATCHES "init")
      list(APPEND HALMD_PAIR_POTENTIALS ${result})
    endif()
    unset(result)
  endforeach(file)
endif()

# define cached variable, enable all external potentials by default
set(HALMD_EXTERNAL_POTENTIALS ${HALMD_EXTERNAL_POTENTIALS} CACHE STRING
    "List of enabled external potentials")
# define cached variable, enable all pair potentials by default
set(HALMD_PAIR_POTENTIALS ${HALMD_PAIR_POTENTIALS} CACHE STRING
    "List of enabled pair potentials")

# define variables HALMD_WITH_*
foreach(potential ${HALMD_EXTERNAL_POTENTIALS})
  list(FIND "${HALMD_EXTERNAL_POTENTIALS}" ${potential} HALMD_WITH_external_${potential})
endforeach()
foreach(potential ${HALMD_PAIR_POTENTIALS})
  list(FIND "${HALMD_PAIR_POTENTIALS}" ${potential} HALMD_WITH_pair_${potential})
endforeach()

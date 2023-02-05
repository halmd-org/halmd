if(${HALMD_WITH_pair_custom})
  add_executable(test_unit_mdsim_potentials_pair_custom
    custom.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_potentials_pair_custom
      halmd_mdsim_gpu_potentials_pair_custom
      halmd_mdsim_gpu
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_potentials_pair_custom
    halmd_mdsim_host_potentials_pair_custom
    halmd_mdsim
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/potentials/pair/custom/host
    test_unit_mdsim_potentials_pair_custom --run_test=custom_host --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/custom/gpu/float
        test_unit_mdsim_potentials_pair_custom --run_test=custom_gpu_float --log_level=test_suite
      )
    endif()
    if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/custom/gpu/dsfloat
        test_unit_mdsim_potentials_pair_custom --run_test=custom_gpu_dsfloat --log_level=test_suite
      )
    endif()
  endif()
endif()

if(${HALMD_WITH_pair_lennard_jones})
  foreach(truncation ${HALMD_PAIR_POTENTIAL_TRUNCATIONS})

    add_executable(test_unit_mdsim_potentials_pair_lennard_jones_${truncation}
      lennard_jones.cpp
    )
    target_compile_definitions(test_unit_mdsim_potentials_pair_lennard_jones_${truncation}
            PRIVATE "-DTRUNCATION_TYPE=${truncation}")
    if(HALMD_WITH_GPU)
      target_link_libraries(test_unit_mdsim_potentials_pair_lennard_jones_${truncation}
        halmd_mdsim_gpu_potentials_pair_lennard_jones
        halmd_mdsim_gpu
        halmd_algorithm_gpu
        halmd_utility_gpu
      )
    endif()
    target_link_libraries(test_unit_mdsim_potentials_pair_lennard_jones_${truncation}
      halmd_mdsim_host_potentials_pair_lennard_jones
      halmd_mdsim
      ${HALMD_TEST_LIBRARIES}
    )
    add_test(unit/mdsim/potentials/pair/lennard_jones/${truncation}/host
      test_unit_mdsim_potentials_pair_lennard_jones_${truncation} --run_test=lennard_jones_host --log_level=test_suite
    )
    if(HALMD_WITH_GPU)
      if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
        halmd_add_gpu_test(unit/mdsim/potentials/pair/lennard_jones/${truncation}/gpu/float
          test_unit_mdsim_potentials_pair_lennard_jones_${truncation} --run_test=lennard_jones_gpu_float --log_level=test_suite
        )
      endif()
      if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
        halmd_add_gpu_test(unit/mdsim/potentials/pair/lennard_jones/${truncation}/gpu/dsfloat
          test_unit_mdsim_potentials_pair_lennard_jones_${truncation} --run_test=lennard_jones_gpu_dsfloat --log_level=test_suite
        )
      endif()
    endif()
  endforeach()
endif()

if(${HALMD_WITH_pair_modified_lennard_jones})
  add_executable(test_unit_mdsim_potentials_pair_modified_lennard_jones
    modified_lennard_jones.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_potentials_pair_modified_lennard_jones
      halmd_mdsim_gpu_potentials_pair_modified_lennard_jones
      halmd_mdsim_gpu
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_potentials_pair_modified_lennard_jones
    halmd_mdsim_host_potentials_pair_modified_lennard_jones
    halmd_mdsim
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/potentials/pair/modified_lennard_jones/host
    test_unit_mdsim_potentials_pair_modified_lennard_jones --run_test=modified_lennard_jones_host --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/modified_lennard_jones/gpu/float
        test_unit_mdsim_potentials_pair_modified_lennard_jones --run_test=modified_lennard_jones_gpu_float --log_level=test_suite
      )
    endif()
    if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/modified_lennard_jones/gpu/dsfloat
        test_unit_mdsim_potentials_pair_modified_lennard_jones --run_test=modified_lennard_jones_gpu_dsfloat --log_level=test_suite
      )
    endif()
  endif()
endif()

if(${HALMD_WITH_pair_power_law})
  add_executable(test_unit_mdsim_potentials_pair_power_law
    power_law.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_potentials_pair_power_law
      halmd_mdsim_gpu_potentials_pair_power_law
      halmd_mdsim_gpu
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_potentials_pair_power_law
    halmd_mdsim_host_potentials_pair_power_law
    halmd_mdsim
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/potentials/pair/power_law/host
    test_unit_mdsim_potentials_pair_power_law --run_test=power_law_host --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/power_law/gpu/float
        test_unit_mdsim_potentials_pair_power_law --run_test=power_law_gpu_float --log_level=test_suite
      )
    endif()
    if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/power_law/gpu/dsfloat
        test_unit_mdsim_potentials_pair_power_law --run_test=power_law_gpu_dsfloat --log_level=test_suite
      )
    endif()
  endif()

  add_executable(test_unit_mdsim_potentials_pair_power_law_hard_core
    power_law_hard_core.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_potentials_pair_power_law_hard_core
      halmd_mdsim_gpu_potentials_pair_power_law
      halmd_mdsim_gpu
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_potentials_pair_power_law_hard_core
    halmd_mdsim_host_potentials_pair_power_law
    halmd_mdsim
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/potentials/pair/power_law_hard_core/host
    test_unit_mdsim_potentials_pair_power_law_hard_core --run_test=power_law_hard_core_host --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/power_law_hard_core/gpu/float
        test_unit_mdsim_potentials_pair_power_law_hard_core --run_test=power_law_hard_core_gpu_float --log_level=test_suite
      )
    endif()
    if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/power_law_hard_core/gpu/dsfloat
        test_unit_mdsim_potentials_pair_power_law_hard_core --run_test=power_law_hard_core_gpu_dsfloat --log_level=test_suite
      )
    endif()
  endif()
endif()

if(${HALMD_WITH_pair_morse})
  add_executable(test_unit_mdsim_potentials_pair_morse
    morse.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_potentials_pair_morse
      halmd_mdsim_gpu_potentials_pair_morse
      halmd_mdsim_gpu
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_potentials_pair_morse
    halmd_mdsim_host_potentials_pair_morse
    halmd_mdsim
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/potentials/pair/morse/host
    test_unit_mdsim_potentials_pair_morse --run_test=morse_host --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/morse/gpu/float
        test_unit_mdsim_potentials_pair_morse --run_test=morse_gpu_float --log_level=test_suite
      )
    endif()
    if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
      halmd_add_gpu_test(unit/mdsim/potentials/pair/morse/gpu/dsfloat
        test_unit_mdsim_potentials_pair_morse --run_test=morse_gpu_dsfloat --log_level=test_suite
      )
    endif()
  endif()
endif()

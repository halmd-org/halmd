foreach(module
    ctype
    error
    file
    id
    )
  add_executable(test_h5xx_${module}
    ${module}.cpp
    )
  target_link_libraries(test_h5xx_${module}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${HDF5_LIBRARY}
    )
  set_target_properties(test_h5xx_${module} PROPERTIES
    LINK_FLAGS "-Wl,--as-needed"
    )
  add_test(h5xx_${module}
    test_h5xx_${module}
    )
endforeach(module)

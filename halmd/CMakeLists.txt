# Generate version header file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

if(NOT HALMD_BACKEND_EXECUTABLES)
  add_executable(halmd
    main.cpp
    options.cpp
    sampler.cpp
    script.cpp
    )
  target_link_libraries(halmd
    -Wl,-whole-archive
    ${HALMD_MODULES}
    -Wl,-no-whole-archive
    ${HALMD_COMMON_LIBRARIES}
    )
  install(TARGETS halmd
    RUNTIME DESTINATION bin
    )
endif(NOT HALMD_BACKEND_EXECUTABLES)

add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(observables)
add_subdirectory(random)
add_subdirectory(utility)

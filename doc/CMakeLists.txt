# Generate documentation in HTML and PDF format using Sphinx.

option(HALMD_GENERATE_DOC "Set to OFF to skip documentation generation even if the documentation tools are found." ON)
if(NOT HALMD_GENERATE_DOC)
  return()
endif(NOT HALMD_GENERATE_DOC)

set(GENERATE_DOC TRUE)

# We use the Sphinx documentation generator to render HTML and manual
# pages from the user and reference documentation in ReST format.
find_package(Sphinx QUIET)
if(NOT SPHINX_FOUND)
  message(WARNING "Unable to find Sphinx documentation generator")
  set(GENERATE_DOC FALSE)
endif()
if(${SPHINX_VERSION_STRING} VERSION_LESS "1.4.0")
  message(WARNING "Version 1.4 or later of Sphinx documentation generator required")
  set(GENERATE_DOC FALSE)
endif()

# Sphinx uses latex and dvipng to render LaTeX formulas.
find_package(LATEX QUIET)
if(NOT LATEX_COMPILER)
  message(WARNING "Unable to find LaTeX compiler")
  set(GENERATE_DOC FALSE)
endif(NOT LATEX_COMPILER)
find_program(DVIPNG_CONVERTER dvipng)
if(NOT DVIPNG_CONVERTER)
  message(WARNING "Unable to find dvipng converter")
  set(GENERATE_DOC FALSE)
endif(NOT DVIPNG_CONVERTER)
mark_as_advanced(DVIPNG_CONVERTER)

# Lua scripts are used to extract inline documentation
find_package(LuaInterp QUIET)
if(NOT LUA_EXECUTABLE)
  message(WARNING "Unable to find Lua interpreter")
  set(GENERATE_DOC FALSE)
endif()

if(NOT GENERATE_DOC)
  message(WARNING "Missing required documentation tools")
  return()
endif()

# Accumulate list of input dependency files for documentation targets,
# to avoid rebuilding documentation when input files remain unchanged.
set(HALMD_DOC_INPUT_FILES)

configure_file(conf.py.in conf.py)
list(APPEND HALMD_DOC_INPUT_FILES ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

# The documentation consists of static ReST documents in the source tree, and
# auto-generated ReST documents in the binary tree. Since Sphinx only supports
# a single tree for sources, copy the static documents to the binary tree.
file(GLOB_RECURSE files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.rst.in")
foreach(input ${files})
  string(REGEX REPLACE "\\.in$" "" output ${input})
  configure_file(${input} ${output} @ONLY)
  list(APPEND HALMD_DOC_INPUT_FILES ${CMAKE_CURRENT_BINARY_DIR}/${output})
endforeach(input)

# Copy image files to the binary tree.
file(GLOB_RECURSE files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "_images/*")
foreach(input ${files})
  if(NOT input MATCHES "(^|/)\\.[^/]+$")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${input}" "${input}" COPYONLY)
    list(APPEND HALMD_DOC_INPUT_FILES ${CMAKE_CURRENT_BINARY_DIR}/${input})
  endif()
endforeach()

# Copy examples for inlining of plot examples into documentation.
file(GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/examples/*)
foreach(input ${files})
  if(NOT input MATCHES "(^|/)\\.[^/]+$")
    configure_file(${CMAKE_SOURCE_DIR}/${input} ${input} @ONLY)
    list(APPEND HALMD_DOC_INPUT_FILES ${CMAKE_CURRENT_BINARY_DIR}/${input})
  endif()
endforeach()

# Copy test scripts for inlining of Lua examples into documentation.
file(GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/test/lua/*)
foreach(input ${files})
  if(NOT input MATCHES "(^|/)\\.[^/]+$")
    configure_file(${CMAKE_SOURCE_DIR}/${input} ${input} @ONLY)
    list(APPEND HALMD_DOC_INPUT_FILES ${CMAKE_CURRENT_BINARY_DIR}/${input})
  endif()
endforeach()

# Extract module documentation from Lua sources
configure_file(modules.lua.in modules.lua)
file(GLOB_RECURSE files RELATIVE "${CMAKE_SOURCE_DIR}/lua" "${CMAKE_SOURCE_DIR}/lua/halmd/*.lua.in")
execute_process(COMMAND
  "${CMAKE_CURRENT_BINARY_DIR}/modules.lua"
    "${CMAKE_SOURCE_DIR}/lua"
    "${CMAKE_CURRENT_BINARY_DIR}/modules"
    ${files}
)
list(APPEND HALMD_DOC_INPUT_FILES ${CMAKE_CURRENT_BINARY_DIR}/modules.lua)
foreach(input ${files})
  string(REGEX REPLACE "\\.in$" "" output ${input})
  list(APPEND HALMD_DOC_INPUT_FILES ${CMAKE_BINARY_DIR}/lua/${output})
endforeach()

# Configure extracted module files (inplace)
#
# Alternatively, we may use ${CMAKE_BINARY_DIR}/lua as input files for modules.lua,
# but there is a dependency problem between lua/CMakeLists.txt and this one.
file(GLOB_RECURSE files "${CMAKE_CURRENT_BINARY_DIR}/modules/*")
foreach(output ${files})
  configure_file("${output}" "${output}.tmp" @ONLY)
  file(RENAME "${output}.tmp" "${output}")
endforeach()

# Define documentation targets
add_custom_command(OUTPUT html
  COMMAND ${SPHINX_EXECUTABLE}
    -q -b html
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/html"
  COMMAND ${CMAKE_COMMAND}
    -E touch html
  DEPENDS ${HALMD_DOC_INPUT_FILES}
  COMMENT "Build HALMD HTML documentation"
)
add_custom_target(halmd_doc_html ALL
  DEPENDS html
)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
  DESTINATION "share/doc/halmd"
)

add_custom_command(OUTPUT man
  COMMAND ${SPHINX_EXECUTABLE}
    -q -b man
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/man"
  COMMAND ${CMAKE_COMMAND}
    -E touch man
  DEPENDS ${HALMD_DOC_INPUT_FILES}
  COMMENT "Build HALMD man page"
)
add_custom_target(halmd_doc_man ALL
  DEPENDS man
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/man/halmd.1"
  DESTINATION "share/man/man1"
)

if(CMAKE_VERBOSE_MAKEFILE OR HALMD_DOC_ONLY)
  set(HALMD_DOC_LATEX_MODE nonstopmode)
else()
  set(HALMD_DOC_LATEX_MODE batchmode)
endif()

add_custom_command(OUTPUT pdf
  COMMAND ${SPHINX_EXECUTABLE}
    -q -b latex
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/pdf"
  COMMAND ${CMAKE_COMMAND}
    -E touch pdf
  DEPENDS ${HALMD_DOC_INPUT_FILES}
  COMMENT "Generate HALMD LaTeX documentation"
)

add_custom_command(OUTPUT pdf/halmd.pdf
  COMMAND make LATEXOPTS=-interaction=${HALMD_DOC_LATEX_MODE}
  DEPENDS pdf
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pdf"
  COMMENT "Build HALMD PDF documentation"
)

add_custom_target(halmd_doc_pdf ALL
  DEPENDS pdf/halmd.pdf
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdf/halmd.pdf"
  DESTINATION "share/doc/halmd"
)

add_custom_target(halmd_doc)
add_dependencies(halmd_doc
  halmd_doc_html
  halmd_doc_man
  halmd_doc_pdf
)

set_directory_properties(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES "html;man;pdf"
)

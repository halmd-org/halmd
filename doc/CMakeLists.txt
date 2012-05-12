# Generate documentation in HTML and PDF format using Sphinx.

set(GENERATE_DOC TRUE)

# We use the Sphinx documentation generator to render HTML and manual
# pages from the user and reference documentation in ReST format.
find_package(Sphinx QUIET)
if(NOT SPHINX_FOUND)
  message(WARNING "Unable to find Sphinx documentation generator")
  set(GENERATE_DOC FALSE)
endif(NOT SPHINX_FOUND)

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

if(NOT GENERATE_DOC)
  message(WARNING "Missing required documentation tools")
  return()
endif()

configure_file(conf.py.in conf.py)

# The documentation consists of static ReST documents in the source tree, and
# auto-generated ReST documents in the binary tree. Since Sphinx only supports
# a single tree for sources, copy the static documents to the binary tree.
file(GLOB_RECURSE files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.rst.in")
foreach(input ${files})
  string(REGEX REPLACE "\\.in$" "" output ${input})
  configure_file(${input} ${output} @ONLY)
endforeach(input)

# Copy examples for inlining of plot examples into documentation.
file(GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/examples/*)
foreach(input ${files})
  configure_file(${CMAKE_SOURCE_DIR}/${input} ${input} @ONLY)
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

add_custom_target(halmd_doc_html ALL
  ${SPHINX_EXECUTABLE}
    -q -b html
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/html"
  COMMENT "Build HALMD HTML documentation"
)

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
  DESTINATION "share/doc/halmd"
)

add_custom_target(halmd_doc_man ALL
  ${SPHINX_EXECUTABLE}
    -q -b man
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/man"
  COMMENT "Build HALMD man page"
)
add_dependencies(halmd_doc_man
  halmd_doc_html
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/man/halmd.1"
  DESTINATION "share/man/man1"
)

add_custom_target(halmd_doc_pdf ALL
  ${SPHINX_EXECUTABLE}
    -q -b latex
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/pdf"
  COMMENT "Build HALMD PDF documentation"
)
add_custom_command(TARGET halmd_doc_pdf POST_BUILD
  COMMAND ${CMAKE_MAKE_PROGRAM} LATEXOPTS=-interaction=nonstopmode
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pdf"
)
add_dependencies(halmd_doc_pdf
  halmd_doc_man
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdf/halmd.pdf"
  DESTINATION "share/doc/halmd"
)

set_directory_properties(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES ".buildinfo;.doctrees;html;man;pdf"
)

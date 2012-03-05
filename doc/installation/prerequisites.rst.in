.. _prerequisites:

Software prerequisites
**********************

.. toctree::
   :maxdepth: 1

   automatic
   manual

The following software packages are required for building HALMD. For an
automated installation procedure, refer to the next section, :ref:`automatic`.
A detailed step-by-step guide for manual installation is given in section
:ref:`manual`.

* `NVIDIA CUDA toolkit <http://www.nvidia.com/object/cuda_get.html>`_ ≥ 3.1

  .. note::

    Versions ≥ 2.2 are supported as well, but as they are incompatible with the
    GCC compiler ≥ 4.4, we recommend to use at least 3.1 with this bug is
    fixed.

  Please refer to the installation instructions shipped with the toolkit.

* `CMake <http://www.cmake.org/>`_ ≥ 2.8.0 with a patch for native CUDA support

  The build process of HALMD depends on CMake, a cross-platform, open-source
  build system.

  .. note::

    This patch adds *native* CUDA source file compilation and linking support
    to CMake and is not to be confused nor compatible with the CUDA module in
    CMake 2.8.

* `Git <http://git-scm.com/>`_ ≥ 1.5.6.2

  The source code of HALMD is managed by Git, a fast and efficient, distributed
  version control system. Git is available for many operating systems and
  their flavours.

* `Boost C++ Libraries <http://www.boost.org/>`_ ≥ 1.43.0

  In addition, the proposed `Boost.Log <http://boost-log.sourceforge.net/>`_
  library is needed.

  The C++ part of HALMD makes deliberate use of libraries in the Boost C++
  collection.

* `Lua interpreter <http://www.lua.org>`_ ≥ 5.1

  A simulation with HALMD is setup and configured by means of the Lua scripting
  language. The fast and lightweight Lua interpreter is embedded in the HALMD
  executable.

* `Luabind library <http://sourceforge.net/projects/luabind>`_ >= 0.9.1

  Luabind is used to create bindings between C++ objects and Lua modules.

* `HDF5 C++ Library <http://www.hdfgroup.org/HDF5/>`_ >= 1.8

  "HDF5 is a data model, library, and file format for storing and managing
  data. It supports an unlimited variety of datatypes, and is designed for
  flexible and efficient I/O and for high volume and complex data."


To optionally generate documentation in HTML and PDF format:

* `Sphinx documentation generator <http://sphinx.pocoo.org/>`_ ≥ 0.6.1

  "Sphinx is a tool that makes it easy to create intelligent and beautiful
  documentation."

* LaTeX including pdflatex and dvipng

* graphviz


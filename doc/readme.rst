.. _installation:

Installation
************

Software prerequisites
======================

These software packages are required for compilation:

* `NVIDIA CUDA toolkit <http://www.nvidia.com/object/cuda_get.html>`_ >= 1.1

  For Tesla C1060 or GeForce GT200-based cards, **CUDA 2.2 or 2.3 is recommended**.

  With CUDA 2.1 we experienced minor performance degradation on a GeForce GTX 280.

  HALMD has proven to run reliably over 10⁸ or more integration steps using
  CUDA driver version 185.18.36 (with CUDA 2.0 or 2.2 runtime) and 190.42 (with
  CUDA 2.3 runtime). Always look for the newest NVIDIA driver version on the
  `NVIDIA Linux driver website <http://www.nvidia.com/object/unix.html>`_, as
  the CUDA driver accompanying a CUDA toolkit is not updated and may contain
  serious bugs, e.g. that cause hanging GPU kernels.

* `CMake <http://www.cmake.org/>`_ >= 2.8.0 with custom CUDA compiler support patch

  The patched CMake version, which adds native CUDA language support, is
  available at ::

    git://git.colberg.org/gpgpu/cmake-cuda

  To obtain the CUDA patch relative to the newestest supported CMake version ::

    git diff cmake..master

  Alternatively, you can use (or self-compile) these Debian packages:

  * `CMakeCUDA packages for Debian squeeze/sid
    <http://colberg.org/debian/pool/main/c/cmake>`_

  * `CMakeCUDA packages for Ubuntu karmic
    <http://colberg.org/ubuntu/pool/main/c/cmake>`_

  .. note::

     This patch adds *native* CUDA source file compilation and linking support
     to CMake and is not to be confused nor compatible with the CUDA module in
     CMake 2.8.

* `Boost C++ Libraries <http://www.boost.org/>`_ >= 1.42.0

  In addition, the proposed `Boost.Log <http://boost-log.sourceforge.net/>`_
  library is needed, which is acquired with ::

    svn co http://boost-log.svn.sourceforge.net/svnroot/boost-log/trunk/boost-log

  To compile the Boost.Log library, copy the directories ``boost/log`` and
  ``libs/log`` to the respective paths in the Boost source tree and
  `compile Boost
  <http://www.boost.org/doc/libs/1_41_0/more/getting_started/unix-variants.html#easy-build-and-install>`_
  following the standard installation procedure.

  Alternatively, you can use (or self-compile) these Debian packages:

  * `Boost packages with Boost.Log for Debian etch
    <http://colberg.org/debian/pool/main/b/boost1.42>`_

  * `Boost packages with Boost.Log for Ubuntu jaunty
    <http://colberg.org/ubuntu/pool/main/b/boost1.42>`_

* `HDF5 C++ Library <http://www.hdfgroup.org/HDF5/>`_ >= 1.6.6

* `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_ (GSL)

* `Git <http://git-scm.com/>`_ >= 1.5.6.2


To optionally generate documentation in HTML and PDF format:

* `Sphinx documentation generator <http://sphinx.pocoo.org/>`_ >= 0.6.1

* LaTeX including pdflatex and dvipng

* graphviz


Compilation
===========

HALMD uses `CMake <http://www.cmake.org/>`_ to generate its make files, which is
similar to the more commonly used Autotools suite recognisable by the
``configure`` script accompanying a software package, but much faster and much
easier to develop with.

With cmake, out-of-tree compilation is preferred, so we generate the compilation
or build tree in a separate directory. This allows one to have multiple builds
of the same software at the same time, for example a release build with
aggressive optimisation and a debug build including debugging symbols. Note that
the build directory may be a subdirectory in the source tree.

Setting up the cmake build tree
-------------------------------

In the cloned HALMD repository, switch to a new build directory::

  mkdir -p build/release && cd build/release

If the third-party packages are installed in standard locations, run ::

  cmake ../..

This will detect all necessary software, and then generate the make files.

Compilation is done using make, which supports parallel builds::

  nice make -j4

Individual backends may be compiled selectively::

  nice make -j4 halmd halmd_gpu_neighbour halmd_host

Note that due to extensive use of C++ templates, a **single compiler process**
may easily consume **more than 500MB of memory**, so be sure to adjust the
number of parallel builds to the available resources.


Updating the build tree
-----------------------

After checking out a new git commit, **switch to the build directory** (e.g.
``build/release``) and run::

  cmake .

This instructs cmake to regenerate the build tree using the configuration from a
previous cmake run. Then compile with ``make`` as usual.


Setting build parameters
------------------------

Parameters may be passed to cmake as environment variables or cache variables.

Environment variables are prepended to the cmake command::

  CXXFLAGS="-fPIC -Wall" cmake ../..

:doc:`cmake/env_vars`

Cache variables are appended using the -D option::

  cmake -DCMAKE_BUILD_TYPE=Release ../..

:doc:`cmake/cache_vars`

The following example demonstrates how to compile separate, dynamically linked
executables for each backend, which are statically linked to all libraries except the
standard C and C++ libraries::

  CXXFLAGS="-fPIC -Wall"
  NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
  cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DHALMD_BACKEND_EXECUTABLES=TRUE \
      -DBoost_USE_STATIC_LIBS=TRUE \
      -DHDF5_USE_STATIC_LIBS=TRUE \
      -DGSL_USE_STATIC_LIBS=TRUE \
      ../..

The options given here correspond to the default values.

Further cmake configuration
---------------------------

Compilation flags may be configured via CMake's text mode interface::

  ccmake .

To finish configuration, hit "c" and "g" to apply and recompile with make.
Alternatively, you may use CMake's graphical interface::

  cmake-gui .

The following switch displays the actual commands invoked by make::

  CMAKE_VERBOSE_MAKEFILE	ON

On some 64-bit systems, cmake may accidently use a 32-bit library instead of its
64-bit counterpart, which results in linker errors. With Mandriva Linux, the
following adjustments are required in ccmake::

  GSL_CBLAS_LIBRARY		/usr/lib64/libgslcblas.so.0
  GSL_LIBRARY			/usr/lib64/libgsl.so.0


An installation prefix may be specified as following::

  CMAKE_INSTALL_PREFIX		/your/home/directory/usr

The compiled program is then installed into this tree by ::

  make install


Testing
=======

HALMD includes a preliminary test suite, which may be started in the build tree by ::

  ctest


Installation
************

License
=======

HALMD is licensed under the GNU General Public License version 3.

Website
=======

The documentation for the current release is available at
http://colberg.org/research/halmd/.


Getting the source code
=======================

HALMD is maintained in a private `Git <http://git-scm.com/>`_ repository,
which requires a git account at the physics faculty of the LMU.
If you are new to Git or version control in general, the `Git tutorial
<http://www.kernel.org/pub/software/scm/git/docs/gittutorial.html>`_
will get you started.
Former Subversion users may also read the `Git SVN Crash Course
<http://git.or.cz/course/svn.html>`_.
For in-depth documentation, see the `Git User's Manual
<http://www.kernel.org/pub/software/scm/git/docs/user-manual.html>`_.

To checkout the main repository::

  git clone username@git.physik.uni-muenchen.de:/pub/scm/granmat/ljgpu

Updates may retrieved within the cloned repository using::

  git pull

Individual releases may be checked out with::

  git checkout halmd-v0.3.6.1

To revert to the main development branch after a release checkout, use::

  git checkout master


Software prerequisites
======================

These software packages are required for compilation:

* `NVIDIA CUDA toolkit <http://www.nvidia.com/object/cuda_get.html>`_ >= 1.1

  For Tesla C1060 or GeForce GT200-based cards, **CUDA 2.0 is strongly
  recommended**.

  With CUDA 2.1 we experienced minor performance degradation on a GeForce GTX
  280, while CUDA 2.2 using the NVIDIA driver 185.18.08 causes hanging GPU
  kernels for long-time simulations running more than a day.

* `CMake <http://www.cmake.org/>`_ >= 2.6.0 with custom CUDA compiler support patch

  The patched CMake version, which adds native CUDA language support, is
  available at::

    git clone username@git.physik.uni-muenchen.de:/pub/scm/granmat/cmake-cuda

* `Boost C++ Libraries <http://www.boost.org/>`_ >= 1.37.0

* `HDF5 C++ Library <http://www.hdfgroup.org/HDF5/>`_ >= 1.6.6

  It is recommended to install a 1.6.x version (instead of 1.8.x), as the
  current version of PyTables is not compatible with HDF5 1.8.x.  Note that
  PyTables is **not** required for HALMD, but for the ``mdplot`` plot package.

* `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_ (GSL)

* `Git <http://git-scm.com/>`_ >= 1.5.6.2

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

If the third-party packages are installed in standard locations, run::

  CXXFLAGS="-fPIC -Wall" \
  CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
  NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
  cmake -DCMAKE_BUILD_TYPE=Release ../..

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

Cache variables are appended using the -D option::

  cmake -DCMAKE_BUILD_TYPE=Release ../..


Useful cmake environment variables
----------------------------------

.. glossary::

   CXXFLAGS
     Compilation flags for C++ compiler.

     Recommended value is ``CXXFLAGS="-fPIC -Wall"``.

   NVCCFLAGS
     Compilation flags for CUDA compiler.

     Recommended value is ``NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c"``.

   CUDA_INSTALL_PREFIX
     Path to CUDA toolkit.

     This variable is not needed if the NVCC compiler is installed in a default
     executable path, e.g. /usr/bin or /usr/local/bin.

   BOOST_LIBRARYDIR
     Path to boost libraries.

   BOOST_INCLUDEDIR
     Path to boost header files.


Useful cmake cache variables
----------------------------

.. glossary::

   CMAKE_BUILD_TYPE
     CMake build type.

     For production builds with -O3 optimisation,
     use ``-DCMAKE_BUILD_TYPE=Release``.

     For debugging builds with -O2 optimisation and debug symbols,
     use ``-DCMAKE_BUILD_TYPE=RelWithDebInfo``.

     To build with CUDA device emulation,
     use ``-DCMAKE_BUILD_TYPE=DeviceEmu``.

   CMAKE_PREFIX_PATH
     Path to third-party libraries, e.g. ``-DCMAKE_PREFIX_PATH=$HOME/usr``.

     This variable is only needed if libraries are installed in non-standard paths.

   Boost_USE_STATIC_LIBS
     Link to Boost libraries statically.

     Recommended value is ``-DBoost_USE_STATIC_LIBS=TRUE``.

   HDF5_USE_STATIC_LIBS
     Link to HDF5 libraries statically.

     Recommended value is ``-DHDF5_USE_STATIC_LIBS=TRUE``.

   GSL_USE_STATIC_LIBS
     Link to GNU Scientific Library statically.

     Recommended value is ``-DGSL_USE_STATIC_LIBS=TRUE``.

   HALMD_BACKEND_EXECUTABLES
     Compile separate, dynamically linked executable for each backend.

     Recommended value is ``DHALMD_BACKEND_EXECUTABLES=TRUE``.

   HALMD_USE_STATIC_LIBS
     Compile separate, statically linked executable for each backend.

     This only compiles the host backends, as the CUDA runtime library requires
     dynamic linking to load the CUDA driver.


   HALMD_VARIANT_CELL_SUMMATION_ORDER
     Use opposite cell summation order (GPU backends).

     Default value is ``TRUE``.

   HALMD_VARIANT_FORCE_DSFUN
     Use double-single precision functions in cell summation (GPU backends).

     Default value is ``TRUE``.

   HALMD_VARIANT_HILBERT_ALT_3D
     Use alternative 3D Hilbert curve vertex rules (GPU backends).

     Default value is ``FALSE``.

   HALMD_VARIANT_HILBERT_ORDER
     Use Hilbert space-filling curve particle ordering (GPU backends).

     Default value is ``TRUE``.

   HALMD_VARIANT_HOST_SINGLE_PRECISION
     Use single-precision math in host implementation (host backends).

     Default value is ``FALSE``.

     This option requires SSE, which is enabled by default on x86_64.

   HALMD_VARIANT_VERLET_DSFUN
     Use double-single precision functions in Verlet integrator (GPU backends).

     Default value is ``TRUE``.


Build parameters example
------------------------

The following example demonstrates how to compile separate, dynamically linked
executables for each backend, and statically link to all libraries except the
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


Further cmake configuration
---------------------------

Compilation flags may be configured via the CMake GUI::

  ccmake .

To finish configuration, hit "c" and "g" to apply and recompile with make.

To display the actual compilation commands, set::

  CMAKE_VERBOSE_MAKEFILE	ON

On some 64-bit systems, cmake may accidently use a 32-bit library instead of its
64-bit counterpart, which results in linker errors. With Mandriva Linux, the
following adjustments are required in ccmake::

  GSL_CBLAS_LIBRARY		/usr/lib64/libgslcblas.so.0
  GSL_LIBRARY			/usr/lib64/libgsl.so.0


An installation prefix may be specified as following::

  CMAKE_INSTALL_PREFIX		/home/Peter.Colberg/usr

The compiled program is then installed into this tree with::

  make install


Testing
=======

HALMD includes a preliminary test suite, which may be started in the build tree with::

  ctest



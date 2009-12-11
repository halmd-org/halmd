Installation
************

Website
=======

The documentation for the current release is available at
http://colberg.org/research/halmd.

Getting the source code
=======================

HALMD is maintained in a `Git <http://git-scm.com/>`_ repository.

If you are new to Git or version control in general, the `Git tutorial
<http://www.kernel.org/pub/software/scm/git/docs/gittutorial.html>`_
will get you started.
Former Subversion users may also read the `Git SVN Crash Course
<http://git.or.cz/course/svn.html>`_.
For in-depth documentation, see the `Git User's Manual
<http://www.kernel.org/pub/software/scm/git/docs/user-manual.html>`_.

Checkout of the main repository ::

  git clone git://git.colberg.org/research/halmd

Updates may be retrieved within the cloned repository using ::

  git pull

A specific version may be checked out with ::

  git checkout halmd-v0.1

To revert to the main development branch after a version checkout, use ::

  git checkout master

License
=======

HALMD is licensed under the GNU General Public License version 3.

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

If the third-party packages are installed in standard locations, run ::

  CXXFLAGS="-fPIC -Wall" \
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

:doc:`env_vars`

Cache variables are appended using the -D option::

  cmake -DCMAKE_BUILD_TYPE=Release ../..

:doc:`cache_vars`

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


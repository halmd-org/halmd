Installation
************

License
=======

ljgpu is licensed under the GNU General Public License version 3.

Website
=======

The documentation for the current release is available at
http://colberg.org/code/ljgpu/.


Getting the source code
=======================

ljgpu is maintained in a private `Git <http://git-scm.com/>`_ repository, which
requires a git account at the physics faculty of the LMU.
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

  git checkout ljgpu-v0.3.6.1

To revert to the main development branch after a release checkout, use::

  git checkout master


Software prerequisites
======================

ljgpu is based on variety of high-quality open source software packages.
If you share the disdain of many scientists of installing third-party software
oneself and are in the fortunate position to freely choose the operating system
on your CUDA machines, you may save yourself a lot of trouble by choosing an
operating system which includes the software in its official package repository,
for example `Debian GNU/Linux <http://www.debian.org/>`_ or `Ubuntu
<http://www.ubuntu.com>`_. The Debian distribution and its derivatives have a
strong focus on universality and therefore come bundled with a plethora of
scientific packages.

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

    git clone username@git.physik.uni-muenchen.de:/pub/scm/granmat/ljgpu

* `Boost C++ Libraries <http://www.boost.org/>`_ >= 1.37.0

* `HDF5 C++ Library <http://www.hdfgroup.org/HDF5/>`_ >= 1.6.6

  It is recommended to install a 1.6.x version (instead of 1.8.x), as the
  current version of PyTables is not compatible with HDF5 1.8.x.  Note that
  PyTables is **not** required for ljgpu, but for the ``mdplot`` plot package.

* `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_ (GSL)

* `Git <http://git-scm.com/>`_ >= 1.5.6.2

* `Sphinx documentation generator <http://sphinx.pocoo.org/>`_ >= 0.6.1

* LaTeX including pdflatex and dvipng


Compilation
===========

ljgpu uses `CMake <http://www.cmake.org/>`_ to generate its make files, which is
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

In the cloned ljgpu repository, switch to a new build directory::

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

  nice make -j4 ljgpu ljgpu_gpu_neighbour ljgpu_host

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


Non-standard third-party library locations
------------------------------------------

If cmake should look for third-party libraries in a custom path, set the
CMAKE_PREFIX_PATH variable::

  CXXFLAGS="-fPIC -Wall" \
  CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
  NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/usr ../..


On non-Debian systems, the Boost header and library locations have to be
passed to cmake::

  CXXFLAGS="-fPIC -Wall" \
  CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
  NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
  BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
  BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
  cmake -DCMAKE_BUILD_TYPE=Release ../..


Static linking
--------------

For static linking to the Boost, HDF5 and GSL libraries, use::

  CXXFLAGS="-fPIC -Wall"
  CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
  NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
  BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
  BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
  cmake -DCMAKE_BUILD_TYPE=Release \
  -DBoost_USE_STATIC_LIBS=TRUE \
  -DHDF5_USE_STATIC_LIBS=TRUE \
  -DGSL_USE_STATIC_LIBS=TRUE \
  ../..

To compile separate, dynamically linked executables for each backend::

  CXXFLAGS="-fPIC -Wall" \
  CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
  NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
  BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
  BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
  cmake -DCMAKE_BUILD_TYPE=Release -Dljgpu_STATIC_BACKEND=TRUE ../..

To compile separate, statically linked executables for each backend::

  CXXFLAGS="-fPIC -Wall" \
  BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
  BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
  cmake -DCMAKE_BUILD_TYPE=Release -Dljgpu_USE_STATIC_LIBS=TRUE ../..

This only compiles the host backends, as the CUDA runtime library
requires dynamic linking to load the CUDA driver.


Debug build
-----------

To configure a build with debugging symbols, switch to a new subdirectory (e.g.
``build/debug``) and run::

  CXXFLAGS="-fPIC -Wall" \
  CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
  NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
  cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../..


Further cmake configuration
---------------------------

Compilation flags may be configured via the CMake GUI::

  ccmake .

To finish configuration, hit "c" and "g" to apply and recompile with make.

The following compilation flags are highly recommended::

  CMAKE_CUDA_FLAGS		--host-compilation=c -Xcompiler -fPIC -Xptxas -v
  CMAKE_CXX_FLAGS		-Wall -fPIC

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

ljgpu includes a preliminary test suite, which may be started in the build tree with::

  ctest



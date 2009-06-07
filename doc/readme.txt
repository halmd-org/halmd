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

To compile, switch to a newly created directory and execute::

  $ CXXFLAGS="-fPIC -Wall" \
	CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
	NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
	cmake -DCMAKE_BUILD_TYPE=Release ~/path/to/ljgpu_source
  $ make

If cmake should look for third-party libraries in a custom path, add::

  -DCMAKE_PREFIX_PATH=$HOME/usr

or ::

  -DCMAKE_INCLUDE_PATH=$HOME/usr/include
  -DCMAKE_LIBRARY_PATH=$HOME/usr/lib

to the arguments passed to cmake and use the following CXXFLAGS::

  CXXFLAGS="-fPIC -Wall -I$HOME/usr/include"


On non-Debian systems, the paths to the Boost headers and libraries
have to be specified, using BOOST_INCLUDEDIR and BOOST_LIBRARYDIR::

  $ CXXFLAGS="-fPIC -Wall" \
	CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
	NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
	BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
	BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
	cmake -DCMAKE_BUILD_TYPE=Release ~/path/to/ljgpu_source
  $ make


For static linking to the Boost, HDF5 and GSL libraries, use::

  $ CXXFLAGS="-fPIC -Wall" \
	CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
	NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
	BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
	BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
	cmake -DCMAKE_BUILD_TYPE=Release \
	-DBoost_USE_STATIC_LIBS=TRUE \
	-DHDF5_USE_STATIC_LIBS=TRUE \
	-DGSL_USE_STATIC_LIBS=TRUE \
	~/path/to/ljgpu_source
  $ make

To compile separate executables for each backend::

  $ CXXFLAGS="-fPIC -Wall" \
	CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
	NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
	BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
	BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
	cmake -DCMAKE_BUILD_TYPE=Release \
	-Dljgpu_STATIC_BACKEND=TRUE \
	~/path/to/ljgpu_source
  $ make

To compile statically linked backend executables::

  $ CXXFLAGS="-fPIC -Wall" \
	BOOST_LIBRARYDIR=/usr/lib64/boost1_37 \
	BOOST_INCLUDEDIR=/usr/include/boost-1_37 \
	cmake -DCMAKE_BUILD_TYPE=Release \
	-Dljgpu_USE_STATIC_LIBS=TRUE \
	~/path/to/ljgpu_source
  $ make

Note that this only compiles the host backends, as the CUDA runtime library
requires dynmamic linking as it dynamically loads the CUDA driver library.

For a build with debugging symbols, use the following::

  $ CXXFLAGS="-fPIC -Wall" \
	CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
	NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ~/path/to/ljgpu_source
  $ make


Compilation flags may be configured via the CMake GUI::

  $ ccmake .

To finish configuration, hit "c" and "g" to apply and recompile with make.

The following compilation flags are highly recommended::

    CMAKE_CUDA_FLAGS		--host-compilation=c -Xcompiler -fPIC -Xptxas -v
    CMAKE_CXX_FLAGS		-Wall -fPIC

To display the actual compilation commands, set::

    CMAKE_VERBOSE_MAKEFILE	ON

On some 64-bit systems (e.g. Mandriva), cmake may accidently use a 32-bit
library instead of its 64-bit counterpart, which results in linker errors.
With Mandriva Linux, the following adjustments are required in ccmake::

    GSL_CBLAS_LIBRARY		/usr/lib64/libgslcblas.so.0
    GSL_LIBRARY			/usr/lib64/libgsl.so.0


An installation prefix may be specified as following::

    CMAKE_INSTALL_PREFIX	/home/Peter.Colberg/usr

The compiled program is then installed into this tree with::

    $ make install


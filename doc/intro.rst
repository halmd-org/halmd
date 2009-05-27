Introduction
************

License
=======

ljgpu is licensed under the GNU General Public License version 3.


Prerequisites
=============

The following software packages are required for compilation.

* NVIDIA CUDA toolkit >= 1.1
* CMake >= 2.6.0 with custom CUDA compiler support patch
* Boost C++ Libraries >= 1.37.0
* HDF5 C++ Library >= 1.6.6
* GNU Scientific Library
* Git >= 1.5.6.2


Installation
============

To compile, switch to a newly created directory and execute::

  $ CXXFLAGS="-fPIC -Wall" \
	CUDA_INSTALL_PREFIX=/path/to/cuda/toolkit \
	NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v --host-compilation=c" \
	cmake -DCMAKE_BUILD_TYPE=Release ~/path/to/ljgpu_source
  $ make

If cmake should look for third-party libraries in a custom path, add::

  -DCMAKE_PREFIX_PATH=$HOME/usr

or

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


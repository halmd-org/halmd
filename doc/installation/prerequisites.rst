.. _prerequisites:

Prerequisites
*************

This section is a step-by-step guide to installing the necessary dependencies
to compile HALMD from source. Be sure to check if your distribution ships with
any of these packages before attempting to compile them yourself. Before
proceeding, be aware of the section :ref:`packages`.

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

* `Git <http://git-scm.com/>`_ >= 1.5.6.2


To optionally generate documentation in HTML and PDF format:

* `Sphinx documentation generator <http://sphinx.pocoo.org/>`_ >= 0.6.1

* LaTeX including pdflatex and dvipng

* graphviz

Overview
========

  * :ref:`prerequisites-cuda` ≥ 3.1
  * :ref:`prerequisites-cmake` ≥ 2.8 with native CUDA support
  * :ref:`prerequisites-boost` ≥ 1.43, supplemented by Boost.Log
  * :ref:`prerequisites-lua` ≥ 5.1
  * :ref:`prerequisites-luabind` ≥ 0.9
  * :ref:`prerequisites-hdf5` ≥ 1.8
  * :ref:`prerequisites-sphinx`

.. tip::

   When installing third-party packages, it is advisable to put them into
   separate directories. If you install software only for yourself, use package
   directories of the form ``~/usr/PKGNAME-PKGVERSION``, for example
   ``~/usr/boost-1.46.1`` or ``~/usr/Sphinx-1.0.4``. If you install software
   system-wide as the root user, use ``/opt/PKGNAME-PKGVERSION``.
   This simple scheme allows you to have multiple versions of a package, or
   remove a package without impacting others.

When initially creating the CMake build tree, include all third-party package
directories in the environment variable ``CMAKE_PREFIX_PATH``.
For example, if Boost, Lua and Luabind are installed in your home directory,
CUDA is installed system-wide, and the HALMD source is in ``~/projects/halmd``,
the initial cmake command might look like this ::

   CMAKE_PREFIX_PATH=~/usr/boost_1_46_1:/opt/cuda-3.1:~/usr/lua-5.1.4:~/usr/luabind-0.9.1 cmake ~/projects/halmd

Instead of setting CMAKE_PREFIX_PATH manually, you would include the package directories in your ~/.bashrc (or your favourite shell's equivalent) ::

   export CMAKE_PREFIX_PATH="${HOME}/usr/boost_1_46_1${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"
   export CMAKE_PREFIX_PATH="/opt/cuda-3.1${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"
   export CMAKE_PREFIX_PATH="${HOME}/usr/lua-5.1.4${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"
   export CMAKE_PREFIX_PATH="${HOME}/usr/luabind-0.9.1${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"

GNU/Linux
=========

.. _prerequisites-cuda:

CUDA
----

.. _prerequisites-cmake:

CMake
-----

The build process of HALMD depends on CMake, a cross-platform, open-source build system.

Get the latest `CMake source package`_, currently `CMake 2.8.4`_.

Get the latest `CMake-CUDA`_ patch, currently `CMake-CUDA 2.8.4-1`_.

To fix a bug with `lib64 paths in CXX-only projects`_, get this `CMake patch`_.

.. _CMake source package: http://cmake.org/cmake/resources/software.html

.. _CMake 2.8.4: http://www.cmake.org/files/v2.8/cmake-2.8.4.tar.gz

.. _CMake-CUDA: http://software.colberg.org/projects/cmake-cuda

.. _CMake-CUDA 2.8.4-1: https://software.colberg.org/attachments/download/20/CMake-CUDA-2.8.4-1-g2ed3c7a.patch

.. _lib64 paths in CXX-only projects: http://public.kitware.com/Bug/view.php?id=10813#c25824

.. _CMake patch: http://public.kitware.com/Bug/file_download.php?file_id=3759&type=bug

Extract the CMake source package, and apply the patches in the CMake source directory with ::

   patch -p1 < ../cmake-cuda-2.8.4.patch
   patch -p1 < ../cmake_add_lib64_prefix_for_any_language.patch

Prepare the CMake build with ::

   ./configure --prefix=$HOME/usr/cmake-cuda-2.8.4

Compile CMake with ::

   make

Install CMake into your packages directory::

   make install

Include CMake in your shell environment, by adding to ~/.bashrc::

   export PATH="${HOME}/usr/cmake-cuda-2.8.4/bin${PATH+:$PATH}"
   export MANPATH="${HOME}/usr/cmake-cuda-2.8.4/man${MANPATH+:$MANPATH}"


.. _prerequisites-boost:

Boost C++ libraries
-------------------

The C++ part of HALMD makes use of a variety of libraries in the Boost C++ collection.

Get the latest `Boost source package`_, currently `Boost 1.46.1`_.

.. _Boost source package: http://www.boost.org/users/download
.. _Boost 1.46.1: http://sourceforge.net/projects/boost/files/boost/1.46.1/boost_1_46_1.tar.bz2

Get the latest `Boost.Log source package`_ from the upstream repository.

.. note::

   The `Boost.Log`_ library is a proposed extension to the Boost C++ libraries.
   As a result of the `formal review of Boost.Log`_ in March 2010, the library has
   been accepted subject to several conditions. It is not shipped yet with
   upstream Boost.

.. _Boost.Log source package: http://boost-log.svn.sourceforge.net/viewvc/boost-log/trunk/boost-log/?view=tar
.. _Boost.Log: http://boost-log.sourceforge.net/
.. _formal review of Boost.Log: http://lists.boost.org/boost-announce/2010/03/0256.php

We will build Boost and Boost.Log in a single step, therefore extract both
source packages and copy the Boost.Log headers and library sources to the
Boost source directory using ::

   cp -r boost-log/boost/log boost_1_46_1/boost/
   cp -r boost-log/libs/log boost_1_46_1/libs/

In the Boost source directory, bootstrap the build with ::

   ./bootstrap.sh

Compile Boost using ::

   ./bjam cxxflags=-fPIC

This compiles both dynamic and static libraries.

.. note:: By default, CMake uses the dynamically linked Boost libraries.

   This is the recommended way of linking to Boost, as static linking of
   the unit test executables significantly increases the size of the build
   tree. If you wish to link statically nevertheless, for example to run a
   program on another machine without your Boost libraries, invoke cmake
   with ``-DBoost_USE_STATIC_LIBS=True`` on the *first* run.

.. warning:: Boost may require more than fifteen minutes to compile.

   You are strongly advised to take a coffee break.

Install the Boost libraries into your packages directory::

   ./bjam cxxflags=-fPIC install --prefix=$HOME/usr/boost_1_46_1

Include Boost in your shell environment, by adding to ~/.bashrc::

   export CMAKE_PREFIX_PATH="${HOME}/usr/boost_1_46_1${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"
   export LD_LIBRARY_PATH="${HOME}/usr/boost_1_46_1/lib${LD_LIBRARY_PATH+:$LD_LIBRARY_PATH}"


.. _prerequisites-lua:

Lua interpreter
---------------

A simulation with HALMD is setup and configured by means of the Lua scripting language. The fast and lightweight Lua interpreter is embedded in the HALMD executable.

Get the latest Lua source package from the `Lua download`_ page, currently `Lua 5.1.4`_.

Get the `Lua 5.1.4-2 patch`_ fixing several bugs.

.. _Lua download: http://www.lua.org/download.html
.. _Lua 5.1.4: http://www.lua.org/ftp/lua-5.1.4.tar.gz
.. _Lua 5.1.4-2 patch: http://www.lua.org/ftp/patch-lua-5.1.4-2

Extract the Lua source package, and apply the patch in the Lua source directory with ::

   cd lua-5.1.4/src
   patch < ../../patch-lua-5.1.4-2

The recommended way of embedding the Lua intepreter in an executable is to link
the Lua library statically, which is the default mode of compilation.

On **32-bit platforms**, compile the Lua library with ::

   make linux

On **64-bit platforms**, include the ``-fPIC`` flag using ::

   make linux CFLAGS="-DLUA_USE_LINUX -fPIC -O2 -Wall"

Install the Lua library into your packages directory::

   make install INSTALL_TOP=~/usr/lua-5.1.4

Include Lua in your shell environment, by adding to ~/.bashrc::

   export CMAKE_PREFIX_PATH="${HOME}/usr/lua-5.1.4${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"
   export PATH="${HOME}/usr/lua-5.1.4/bin${PATH+:$PATH}"
   export MANPATH="${HOME}/usr/lua-5.1.4/man${MANPATH+:$MANPATH}"


.. _prerequisites-luabind:

Luabind library
---------------

Luabind is used to create bindings between C++ objects and Lua modules.

Get the latest `Luabind source package`_, currently `Luabind 0.9.1`_.

.. _Luabind source package: http://sourceforge.net/projects/luabind/files/luabind
.. _Luabind 0.9.1: http://sourceforge.net/projects/luabind/files/luabind/0.9.1/luabind-0.9.1.tar.gz

.. note::

   Luabind is based on the Boost C++ libraries and uses boost-jam as its
   build tool. After bootstrapping Boost following the instructions above, the
   bjam executable is found in the top-level source directory, for example
   ``/tmp/boost_1_46_1/bjam``. This directory also has to be passed to bjam
   during Luabind build using the environment variable ``BOOST_ROOT``.

Compile statically linked release and debug variants of the Luabind library with ::

   BOOST_ROOT=/tmp/boost_1_46_1 LUA_PATH=~/usr/lua-5.1.4 /tmp/boost_1_46_1/bjam cxxflags=-fPIC link=static variant=release variant=debug

Install the Luabind library into your packages directory::

   BOOST_ROOT=/tmp/boost_1_46_1 LUA_PATH=~/usr/lua-5.1.4 /tmp/boost_1_46_1/bjam cxxflags=-fPIC link=static variant=release variant=debug install --prefix=$HOME/usr/luabind-0.9.1

(Note that bjam does not replace ~ with your home directory, use ``$HOME`` instead.)

Include Luabind in your shell environment, by adding to ~/.bashrc::

   export CMAKE_PREFIX_PATH="${HOME}/usr/luabind-0.9.1${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"


.. _prerequisites-hdf5:

HDF5 library
------------

"HDF5 is a data model, library, and file format for storing and managing data.
It supports an unlimited variety of datatypes, and is designed for flexible and
efficient I/O and for high volume and complex data."

Get the latest `HDF5 source package`_, currently `HDF5 1.8.6`_.

.. _HDF5 source package: http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain
.. _HDF5 1.8.6: http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.6.tar.gz

Prepare a statically linked build of the HDF5 C and C++ library with ::

   CFLAGS=-fPIC CXXFLAGS=-fPIC ./configure --enable-cxx --enable-static --disable-shared --prefix=$HOME/usr/hdf5-1.8.6

.. note:: Compiling HDF5 with C++ support disables multi-threading.

Compile HDF5 using ::

   make

Install the HDF5 libraries into your packages directory::

   make install

Include HDF5 in your shell environment, by adding to ~/.bashrc::

   export PATH="${HOME}/usr/hdf5-1.8.6/bin${PATH+:$PATH}"
   export CMAKE_PREFIX_PATH="${HOME}/usr/hdf5-1.8.6${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"


.. _prerequisites-sphinx:

Sphinx documentation generator
------------------------------

"Sphinx is a tool that makes it easy to create intelligent and beautiful
documentation."

Get the latest `Sphinx source package`_, currently `Sphinx 1.0.7`_.

.. _Sphinx source package: http://pypi.python.org/pypi/Sphinx
.. _Sphinx 1.0.7: http://pypi.python.org/packages/source/S/Sphinx/Sphinx-1.0.7.tar.gz

Query your Python version ::

   python -V

Create a package directory for Sphinx using the Python major and minor version ::

   mkdir -p $HOME/usr/Sphinx-1.0.7/lib/python2.5/site-packages

Add the package directory to the PYTHON_PATH environment variable ::

   export PYTHONPATH="${HOME}/usr/Sphinx-1.0.7/lib/python2.5/site-packages${PYTHONPATH+:$PYTHONPATH}"

Install Sphinx into your packages directory ::

   python setup.py install --prefix=$HOME/usr/Sphinx-1.0.7

Include Sphinx in your shell environment, by adding to ~/.bashrc::

   export PATH="${HOME}/usr/Sphinx-1.0.7/bin${PATH+:$PATH}"
   export PYTHONPATH="${HOME}/usr/Sphinx-1.0.7/lib/python2.5/site-packages${PYTHONPATH+:$PYTHONPATH}"

AIX
===

Boost
-----

Compile and install Boost using ::

   ./bjam --toolset=vacpp address-model=64 cxxflags=-qrtti=all install --prefix=$HOME/usr/powerpc-ibm-aix5.3/boost_1_46_1

Lua
---

Compile the Lua library ::

   make aix


Luabind
-------

Compile and install the Luabind library ::

   BOOST_ROOT=/tmp/boost_1_46_1 LUA_PATH=~/usr/powerpc-ibm-aix5.3/lua-5.1.4 /tmp/boost_1_46_1/bjam --toolset=vacpp address-model=64 cxxflags=-qrtti=all link=static variant=release variant=debug install --prefix=$HOME/usr/powerpc-ibm-aix5.3/luabind-0.9.1


HDF5
----

Prepare a statically linked build of the HDF5 C and C++ library with ::

   CC=xlc_r CXX=xlC_r CXXFLAGS=-qrtti=all ./configure --enable-cxx --enable-static --disable-shared --prefix=$HOME/usr/powerpc-ibm-aix5.3/hdf5-1.8.6


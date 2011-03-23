Prerequisites
=============

This section is a step-by-step guide to installing the necessary dependencies to
compile HALMD from source. Be sure to check if your distribution ships with any
of these packages before attempting to compile them yourself.

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


CMake
-----

Get the latest `CMake source package`_, currently `CMake 2.8.4`_.

Get the latest `CMake-CUDA`_ patch, currently `CMake-CUDA 2.8.4`_.

To fix a bug with `lib64 paths in CXX-only projects`_, get this `CMake patch`_.

.. _CMake source package: http://cmake.org/cmake/resources/software.html

.. _CMake 2.8.4: http://www.cmake.org/files/v2.8/cmake-2.8.4.tar.gz

.. _CMake-CUDA: http://software.colberg.org/projects/cmake-cuda

.. _CMake-CUDA 2.8.4: http://software.colberg.org/attachments/download/17/CMake-CUDA-2.8.4-0-g7f92dd3.patch

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


Boost
-----

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

   ./bjam

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

   ./bjam install --prefix=$HOME/usr/boost_1_46_1

Include Boost in your shell environment, by adding to ~/.bashrc::

   export CMAKE_PREFIX_PATH="${HOME}/usr/boost_1_46_1${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"
   export LD_LIBRARY_PATH="${HOME}/usr/boost_1_46_1/lib${LD_LIBRARY_PATH+:$LD_LIBRARY_PATH}"


Lua
---

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

   make linux CFLAGS='-fPIC -O2 -Wall $(MYCFLAGS)'

(Note the single quotes to prevent the shell from swallowing $.)

Install the Lua library into your packages directory::

   make install INSTALL_TOP=~/usr/lua-5.1.4

Include Lua in your shell environment, by adding to ~/.bashrc::

   export CMAKE_PREFIX_PATH="${HOME}/usr/lua-5.1.4${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"
   export PATH="${HOME}/usr/lua-5.1.4/bin${PATH+:$PATH}"
   export MANPATH="${HOME}/usr/lua-5.1.4/man${MANPATH+:$MANPATH}"


Luabind
-------

Get the latest `Luabind source package`_, currently `Luabind 0.9.1`_.

.. _Luabind source package: http://sourceforge.net/projects/luabind/files/luabind
.. _Luabind 0.9.1: http://sourceforge.net/projects/luabind/files/luabind/0.9.1/luabind-0.9.1.tar.gz

.. note::

   Luabind is based on the Boost C++ libraries and uses boost-jam as its
   build tool. After bootstrapping Boost following the instructions above, the
   bjam executable is found in the top-level source directory, for example
   ``/tmp/boost_1_46_1/bjam``. This directory also has to be passed to bjam
   during Luabind build using the environment variable ``BOOST_ROOT``.

Compile a statically linked release build of the Luabind library with ::

   BOOST_ROOT=/tmp/boost_1_46_1 LUA_PATH=~/usr/lua-5.1.4 /tmp/boost_1_46_1/bjam link=static variant=release

Install the Luabind library into your packages directory::

   BOOST_ROOT=/tmp/boost_1_46_1 LUA_PATH=~/usr/lua-5.1.4 /tmp/boost_1_46_1/bjam link=static variant=release install --prefix=$HOME/usr/luabind-0.9.1

(Note that bjam does not replace ~ with your home directory, use ``$HOME`` instead.)

Include Luabind in your shell environment, by adding to ~/.bashrc::

   export CMAKE_PREFIX_PATH="${HOME}/usr/luabind-0.9.1${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"


HDF5
----

Get the latest `HDF5 source package`_, currently `HDF5 1.8.6`_.

.. _HDF5 source package: http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain
.. _HDF5 1.8.6: http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.6.tar.gz

Prepare a statically linked build of the HDF5 C and C++ library with ::

   ./configure --enable-cxx --enable-static --disable-shared --prefix=$HOME/usr/hdf5-1.8.6

.. note:: Compiling HDF5 with C++ support disables multi-threading.

Compile HDF5 using ::

   make

Install the HDF5 libraries into your packages directory::

   make install

Include HDF5 in your shell environment, by adding to ~/.bashrc::

   export PATH="${HOME}/usr/hdf5-1.8.6/bin${PATH+:$PATH}"
   export CMAKE_PREFIX_PATH="${HOME}/usr/hdf5-1.8.6${CMAKE_PREFIX_PATH+:$CMAKE_PREFIX_PATH}"


Sphinx
------

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


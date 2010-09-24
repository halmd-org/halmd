Prerequisites
=============

This section is a step-by-step guide to installing the necessary dependencies to
compile HALMD from source. Be sure to check if your distribution ships with any
of these packages before attempting to compile them yourself.

.. tip::

   When installing third-party packages, it is advisable to put them into
   separate directories. If you install software only for yourself, use package
   directories of the form ``~/usr/PKGNAME-PKGVERSION``, for example
   ``~/usr/boost-1.44.0`` or ``~/usr/Sphinx-1.0.4``. If you install software
   system-wide as the root user, use ``/opt/PKGNAME-PKGVERSION``.
   This simple scheme allows you to have multiple versions of a package, or
   remove a package without impacting others.

When initially creating the CMake build tree, include all third-party package
directories in the CMake variable ``CMAKE_PREFIX_PATH``.
For example, if Boost, Lua and Luabind are installed in your home directory,
CUDA is installed system-wide, and the HALMD source is in ``~/projects/halmd``,
the initial cmake command might look like this ::

   cmake -DCMAKE_PREFIX_PATH='~/usr/boost_1_44_0;/opt/cuda-3.1;~/usr/lua-5.1.4;~/usr/luabind-0.9' ~/projects/halmd

(Note the single quotes to prevent the shell from swallowing semicolons.)


CMake
-----

Get the latest `CMake with native CUDA support`_ ::

   git clone git://git.colberg.org/gpgpu/cmake-cuda.git

.. note::

   If you are trapped behind a restrictive firewall, use ::

      git clone http://git.colberg.org/gpgpu/cmake-cuda.git

.. _CMake with native CUDA support: https://software.colberg.org/projects/cmake-cuda

Prepare the CMake build with ::

   ./configure --prefix=$HOME/usr/cmake-cuda-2.8.2

Compile CMake with ::

   make

Install CMake into your packages directory::

   make install


Boost
-----

Get the latest `Boost source package`_, currently `Boost 1.44.0`_.

.. _Boost source package: http://www.boost.org/users/download
.. _Boost 1.44.0: http://sourceforge.net/projects/boost/files/boost/1.44.0/boost_1_44_0.tar.bz2

Get the latest `Boost.Log source package`_, currently `Boost.Log 1.0`_.

.. note::

   The `Boost.Log`_ library is a proposed extension to the Boost C++ libraries.
   As a result of the `formal review of Boost.Log`_ in March 2010, the library has
   been accepted subject to several conditions. It is not shipped yet with
   upstream Boost.

.. _Boost.Log source package: http://sourceforge.net/projects/boost-log/files
.. _Boost.Log 1.0: http://sourceforge.net/projects/boost-log/files/boost-log-1.0.zip
.. _Boost.Log: http://boost-log.sourceforge.net/
.. _formal review of Boost.Log: http://lists.boost.org/boost-announce/2010/03/0256.php

We will build Boost and Boost.Log in a single step, therefore extract both
source packages and copy the Boost.Log headers and library sources to the
Boost source directory using ::

   cp -r /tmp/boost-log-1.0/boost/log /tmp/boost_1_44_0/boost/
   cp -r /tmp/boost-log-1.0/libs/log /tmp/boost_1_44_0/libs/

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

   ./bjam install --prefix=$HOME/usr/boost_1_44_0


Lua
---

Get the latest Lua source package from the `Lua download`_ page, currently `Lua 5.1.4`_.

.. _Lua download: http://www.lua.org/download.html
.. _Lua 5.1.4: http://www.lua.org/ftp/lua-5.1.4.tar.gz

The recommended way of embedding the Lua intepreter in an executable is to link
the Lua library statically, which is the default mode of compilation.

On **32-bit platforms**, compile the Lua library with ::

   make linux

On **64-bit platforms**, include the ``-fPIC`` flag using ::

   make linux CFLAGS='-fPIC -O2 -Wall $(MYCFLAGS)'

(Note the single quotes to prevent the shell from swallowing $.)

Install the Lua library into your packages directory::

   make install INSTALL_TOP=~/usr/lua-5.1.4


Luabind
-------

Get the latest `Luabind source package`_, currently `Luabind 0.9`_.

.. _Luabind source package: http://sourceforge.net/projects/luabind/files/luabind
.. _Luabind 0.9: http://sourceforge.net/projects/luabind/files/luabind/0.9/luabind-0.9.tar.gz

.. note::

   Luabind is based on the Boost C++ libraries and uses boost-jam as its
   build tool. After bootstrapping Boost following the instructions above, the
   bjam executable is found in the top-level source directory, for example
   ``/tmp/boost_1_44_0/bjam``. This directory also has to be passed to bjam
   during Luabind build using the environment variable ``BOOST_ROOT``.

Compile a statically linked release build of the Luabind library with ::

   BOOST_ROOT=/tmp/boost_1_44_0 LUA_PATH=~/usr/lua-5.1.4 /tmp/boost_1_44_0/bjam link=static variant=release

Install the Luabind library into your packages directory::

   BOOST_ROOT=/tmp/boost_1_44_0 LUA_PATH=~/usr/lua-5.1.4 /tmp/boost_1_44_0/bjam link=static variant=release install --prefix=$HOME/usr/luabind-0.9

(Note that bjam does not replace ~ with your home directory, use ``$HOME`` instead.)


HDF5
----

Get the latest `HDF5 source package`_, currently `HDF5 1.8.5 patch 1`_.

.. _HDF5 source package: http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain
.. _HDF5 1.8.5 patch 1: http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.5-patch1.tar.gz

Prepare a statically linked build of the HDF5 C and C++ library with ::

   ./configure --enable-cxx --enable-static --disable-shared --prefix=$HOME/usr/hdf5-1.8.5-patch1

.. note:: Compiling HDF5 with C++ support disables multi-threading.

Compile HDF5 using ::

   make

Install the HDF5 libraries into your packages directory::

   make install


GNU Scientific Library
----------------------

NVIDIA CUDA toolkit
-------------------

Sphinx
------


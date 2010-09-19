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

Boost
-----

Lua
---

Luabind
-------

HDF5
----

GNU Scientific Library
----------------------

NVIDIA CUDA toolkit
-------------------

Sphinx
------


Installing prerequisites with packages.mk
=========================================

This guide describes an easy way of installing all packages required for HALMD,
so you, a novice HALMD user, may continue right away to installing HALMD itself.
To find out more about the packages used by HALMD, see :ref:`prerequisites`.

Quick Start Guide
-----------------

Create and change to a new directory (preferably on a local disk)::

  mkdir /tmp/halmd_prerequisites && cd /tmp/halmd_prerequisites

Download, compile (with 4 processes) and install required packages to ``~/opt``::

  nice make -f .../halmd/examples/packages.mk CONCURRENCY_LEVEL=4 install

Add packages to shell environment::

  make -f .../halmd/examples/packages.mk env >> ~/.bashrc

That was easy!

A more thorough look at packages.mk
-----------------------------------

The packages makefile provides many more rules than just ``install``::

  make -f .../halmd/examples/packages.mk TAB TAB

::

  build              clean-luabind      env-boost          fetch-boost
  build-boost        configure-boost    env-cmake          fetch-cmake
  build-cmake        configure-cmake    env-hdf5           fetch-hdf5
  build-hdf5         configure-hdf5     env-lua            fetch-lua
  build-lua          distclean          env-luabind        fetch-luabind
  build-luabind      distclean-boost    extract-boost      install
  clean              distclean-cmake    extract-cmake      install-boost
  clean-boost        distclean-hdf5     extract-hdf5       install-cmake
  clean-cmake        distclean-lua      extract-lua        install-hdf5
  clean-hdf5         distclean-luabind  extract-luabind    install-lua
  clean-lua          env                fetch              install-luabind

You may choose to install only selected dependencies::

  make -f .../halmd/examples/packages.mk install-boost install-cmake

To compile and install to a path other than ``~/opt``::

  make -f .../halmd/examples/packages.mk install PREFIX=~/pkg/debian6.0-x86_64

If you wish to first download all packages::

  make -f .../halmd/examples/packages.mk fetch

To remove all package build directories::

  make -f .../halmd/examples/packages.mk clean

Also remove downloaded tarballs and patches::

  make -f .../halmd/examples/packages.mk distclean

To compile as a non-root user and install as root::

  make -f .../halmd/examples/packages.mk PREFIX=/opt
  sudo make -f .../halmd/examples/packages.mk install PREFIX=/opt


Troubleshooting
---------------

There are some requirements to ensure a smooth run of packages.mk:

    - a recent C++ compiler (e.g., GCC)
    - some standard tools: wget, tar, gzip, rm, cp, touch, patch
    - Boost requires a development package of the bzip2 library (files ``bzlib.h`` and ``libbz2.so``)

The bzip2 library is necessary for Boost.IOStreams. As HAL's MD package does
not make use of this library, you may opt to compile Boost without bzip2
support by prepending ``NO_BZIP2=1`` to the make command.

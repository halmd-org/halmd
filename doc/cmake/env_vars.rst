Useful environment variables for CMake
--------------------------------------

.. glossary::

   CXXFLAGS
     Compilation flags for C++ compiler.

     Recommended value is ``CXXFLAGS="-fPIC -Wall"``.

   NVCCFLAGS
     Compilation flags for CUDA compiler.

     Recommended value is ``NVCCFLAGS="-Xcompiler -fPIC -Xptxas -v"``.

   CUDA_INSTALL_PREFIX
     Path to CUDA toolkit.

     This variable is not needed if the NVCC compiler is installed in a default
     executable path, e.g. /usr/bin or /usr/local/bin, or in one of these
     standard locations: /usr/local/cuda, /usr/lib/cuda, /usr/shared/cuda,
     /opt/cuda

   BOOST_LIBRARYDIR
     Path to boost libraries.

   BOOST_INCLUDEDIR
     Path to boost header files.


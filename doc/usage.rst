Usage
*****

Prerequisites
=============

To allow the shell to find the program and its required libraries, set these
environment variables in the shell initialization file (e.g. ~/.bashrc)::

    export PATH="$HOME/usr/bin${PATH+:$PATH}"
    export LD_LIBRARY_PATH="$HOME/usr/lib${LD_LIBRARY_PATH+:$LD_LIBRARY_PATH}"


Getting started
===============

ljgpu has three ways of accepting program parameters:

* pass directly as command line parameters
* read from parameter input file [INI format]
* read from HDF5 data file, optionally resuming from a prior trajectory

Options are described in the command line help::

    $ ljgpu --help


Multi-GPU machines
==================

On a machine with multiple GPUs, it is advisable to use the nvlock
program to exclusively assign a GPU to each process::

    $ nvlock ljgpu [...]

This is only necessary for CUDA <= 2.1.

For optimal exploitation of a multi-GPU machine, a job scheduler is highly
recommended, e.g. SLURM — Simple Linux Utility for Resource Management.
https://computing.llnl.gov/linux/slurm/


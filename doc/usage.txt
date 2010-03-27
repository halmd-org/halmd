Usage
*****

Getting started
===============

Program parameters
------------------

HALMD has three ways of accepting program parameters:

* pass directly as command line parameters
* read from parameter input file [INI format]
* read from HDF5 data file, optionally resuming from a prior trajectory

Options are described in the command line help::

  halmd --help


Simulation backends
-------------------

HALMD implements various simulation backends using dynamically loaded plugins.

.. glossary::

   host backend
     A single-threaded soft-sphere molecular dynamics simulation.

     This is a straight-forward reference implementation using cell and
     neighbour lists.

     D. C. Rapaport,
     *The Art of Molecular Dynamics Simulation*,
     Cambridge University Press, 2004

   gpu_square backend
     A naive soft-sphere molecular dynamics simulation for the GPU.

     This is a straight-forward implementation without cell or neighbour lists.
     Simulation time therefore **scales quadratically** with the system size.

     J. A. van Meel, A. Arnold, D. Frenkel, S.F. Portegies Zwart and R. G. Belleman,
     *Harvesting graphics power for MD simulations*,
     `Simulation, Taylor & Francis, 2008, 34, 259-266 <http://dx.doi.org/10.1080/08927020701744295>`_

   gpu_neighbour backend
     An optimised soft-sphere molecular dynamics simulation for the GPU.

     This implementation for the GPU uses fixed-size cell and neighbour lists, which yields
     a speedup over the CPU of order 100.

     P. H. Colberg and F. Höfling,
     *Accelerating glassy dynamics using graphics processing units*,
     `arXiv:0912.3824 [cond-mat.soft] <http://arxiv.org/abs/0912.3824>`_

     J. A. Anderson, C. D. Lorenz and A. Travesset,
     *General purpose molecular dynamics simulations fully implemented on graphics processing units*,
     `Journal of Computational Physics, 2008, 227, 5342-5359 <http://dx.doi.org/10.1016/j.jcp.2008.01.047>`_

     J. A. van Meel, A. Arnold, D. Frenkel, S.F. Portegies Zwart and R. G. Belleman,
     *Harvesting graphics power for MD simulations*,
     `Simulation, Taylor & Francis, 2008, 34, 259-266 <http://dx.doi.org/10.1080/08927020701744295>`_

     H. Sagan,
     *A three-dimensional Hilbert curve*,
     `International Journal of Mathematical Education in Science and Technology, Taylor & Francis, 1993, 24, 541-545 <http://dx.doi.org/10.1080/0020739930240405>`_


Multi-GPU machines
==================

To distribute multiple HALMD processes among CUDA devices in a single machine,
the CUDA devices have to be locked exclusively by the respective process.
HALMD will then choose the first available CUDA device, or an available device
in the subset given by the ``--device`` option.

nvidia-smi tool
---------------

If your NVIDIA driver version comes with the nvidia-smi tool, set all CUDA
devices to *compute exclusive mode* to restrict use to one process per device::

  sudo nvidia-smi --gpu=0 --compute-mode-rules=1
  sudo nvidia-smi --gpu=1 --compute-mode-rules=1

.. warning::

  Compute exclusive mode seems to work reliably only with NVIDIA Tesla devices.
  Although NVIDIA GeForce cards may be set to compute exclusive mode as well,
  doing so might occasionally cause a system crash.

nvlock tool
-----------

If your NVIDIA driver version does not support the nvidia-smi tool, or if you
wish not to set the devices to compute exclusive mode, the ``nvlock`` tool
may be used to exclusively assign a GPU to each process::

  nvlock halmd [...]

You may also directly use the preload library::

  LD_PRELOAD=libnvlock.so halmd [...]

nvlock is available at ::

  git://git.colberg.org/gpgpu/nvlock

or using the "dumb" HTTP protocol ::

  http://git.colberg.org/gpgpu/nvlock

and is compiled with ::

  cmake .
  make


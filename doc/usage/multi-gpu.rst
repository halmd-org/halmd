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

